from __future__ import annotations

import copy
import re
import shutil
from pathlib import Path

from lxml import etree
from PIL import Image, ImageDraw, ImageFont
from pypdf import PageObject, PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas


SUPPLEMENTARY_FIGURE_GROUPS = {
    1: (1, 2),
    2: (4,),
    3: (5, 6),
    4: (7, 9),
    5: (8,),
    6: (10,),
    7: (12,),
    8: (15,),
}


def _font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for candidate in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf",
    ):
        if Path(candidate).is_file():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def _compose_png(source_paths: list[Path], target: Path) -> None:
    images = []
    for source in source_paths:
        with Image.open(source) as image:
            rgba = image.convert("RGBA")
            background = Image.new("RGBA", rgba.size, "white")
            background.alpha_composite(rgba)
            images.append(background.convert("RGB"))

    target_width = max(image.width for image in images)
    label_height = 100 if len(images) > 1 else 0
    gap = 70 if len(images) > 1 else 0
    scaled: list[Image.Image] = []
    for image in images:
        if image.width != target_width:
            height = round(image.height * target_width / image.width)
            image = image.resize((target_width, height), Image.Resampling.LANCZOS)
        scaled.append(image)

    total_height = sum(image.height + label_height for image in scaled)
    total_height += gap * max(0, len(scaled) - 1)
    canvas_image = Image.new("RGB", (target_width, total_height), "white")
    draw = ImageDraw.Draw(canvas_image)
    panel_font = _font(52)
    y = 0
    for index, image in enumerate(scaled):
        if label_height:
            draw.text((38, y + 22), chr(65 + index), fill="#171C24", font=panel_font)
            y += label_height
        canvas_image.paste(image, (0, y))
        y += image.height + gap
    target.parent.mkdir(parents=True, exist_ok=True)
    canvas_image.save(target, dpi=(400, 400), optimize=True)


def _compose_pdf(source_paths: list[Path], target: Path) -> None:
    pages = [PdfReader(str(source)).pages[0] for source in source_paths]
    if len(pages) == 1:
        shutil.copy2(source_paths[0], target)
        return

    target_width = 720.0
    margin = 24.0
    label_height = 24.0
    gap = 20.0
    rendered: list[tuple[PageObject, float, float]] = []
    for page in pages:
        source_width = float(page.mediabox.width)
        source_height = float(page.mediabox.height)
        scale = (target_width - 2 * margin) / source_width
        rendered.append((page, scale, source_height * scale))
    target_height = (
        2 * margin
        + sum(label_height + height for _, _, height in rendered)
        + gap * (len(rendered) - 1)
    )
    output_page = PageObject.create_blank_page(
        width=target_width,
        height=target_height,
    )

    label_positions: list[tuple[str, float]] = []
    cursor_top = target_height - margin
    for index, (page, scale, height) in enumerate(rendered):
        label_y = cursor_top - 16
        label_positions.append((chr(65 + index), label_y))
        cursor_top -= label_height
        y = cursor_top - height
        transform = Transformation().scale(scale).translate(margin, y)
        output_page.merge_transformed_page(page, transform)
        cursor_top = y - gap

    overlay_path = target.with_suffix(".labels.pdf")
    overlay = canvas.Canvas(str(overlay_path), pagesize=(target_width, target_height))
    overlay.setFont("Helvetica-Bold", 15)
    overlay.setFillColorRGB(0.09, 0.11, 0.14)
    for label, y in label_positions:
        overlay.drawString(margin, y, label)
    overlay.save()
    output_page.merge_page(PdfReader(str(overlay_path)).pages[0])
    writer = PdfWriter()
    writer.add_page(output_page)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as handle:
        writer.write(handle)
    overlay_path.unlink(missing_ok=True)


def _svg_dimensions(root: etree._Element) -> tuple[float, float]:
    view_box = root.get("viewBox")
    if view_box:
        values = [float(value) for value in re.split(r"[,\s]+", view_box.strip())]
        if len(values) == 4:
            return values[2], values[3]

    def numeric(value: str | None) -> float:
        if not value:
            return 0.0
        match = re.search(r"[-+]?\d*\.?\d+", value)
        return float(match.group()) if match else 0.0

    width = numeric(root.get("width"))
    height = numeric(root.get("height"))
    if width <= 0 or height <= 0:
        raise ValueError("SVG lacks a usable viewBox or width/height")
    return width, height


def _prefix_svg_ids(element: etree._Element, prefix: str) -> None:
    mapping: dict[str, str] = {}
    for node in element.iter():
        identifier = node.get("id")
        if identifier:
            mapping[identifier] = f"{prefix}{identifier}"
            node.set("id", mapping[identifier])
    if not mapping:
        return
    for node in element.iter():
        for attribute, value in list(node.attrib.items()):
            updated = value
            for old, new in mapping.items():
                updated = updated.replace(f"url(#{old})", f"url(#{new})")
                if updated == f"#{old}":
                    updated = f"#{new}"
            if updated != value:
                node.set(attribute, updated)
        if isinstance(node.tag, str) and node.tag.endswith("style") and node.text:
            for old, new in mapping.items():
                node.text = node.text.replace(f"url(#{old})", f"url(#{new})")


def _compose_svg(source_paths: list[Path], target: Path) -> None:
    if len(source_paths) == 1:
        shutil.copy2(source_paths[0], target)
        return

    parsed: list[tuple[etree._Element, float, float]] = []
    for index, source in enumerate(source_paths):
        root = etree.parse(str(source)).getroot()
        width, height = _svg_dimensions(root)
        imported = copy.deepcopy(root)
        _prefix_svg_ids(imported, f"panel{index}_")
        parsed.append((imported, width, height))

    target_width = 3600.0
    margin = 80.0
    label_height = 100.0
    gap = 80.0
    rendered_heights = [
        height * (target_width - 2 * margin) / width
        for _, width, height in parsed
    ]
    target_height = (
        2 * margin
        + sum(label_height + height for height in rendered_heights)
        + gap * (len(parsed) - 1)
    )
    svg_ns = "http://www.w3.org/2000/svg"
    outer = etree.Element(
        f"{{{svg_ns}}}svg",
        nsmap={None: svg_ns},
        width=f"{target_width:g}",
        height=f"{target_height:g}",
        viewBox=f"0 0 {target_width:g} {target_height:g}",
    )
    etree.SubElement(
        outer,
        f"{{{svg_ns}}}rect",
        x="0",
        y="0",
        width=f"{target_width:g}",
        height=f"{target_height:g}",
        fill="#ffffff",
    )

    cursor = margin
    for index, ((imported, width, _), rendered_height) in enumerate(
        zip(parsed, rendered_heights)
    ):
        label = etree.SubElement(
            outer,
            f"{{{svg_ns}}}text",
            x=f"{margin:g}",
            y=f"{cursor + 58:g}",
            fill="#171c24",
            style="font-family:Arial,sans-serif;font-size:52px;font-weight:700",
        )
        label.text = chr(65 + index)
        cursor += label_height
        scale = (target_width - 2 * margin) / width
        group = etree.SubElement(
            outer,
            f"{{{svg_ns}}}g",
            transform=f"translate({margin:g},{cursor:g}) scale({scale:g})",
        )
        for child in imported:
            group.append(copy.deepcopy(child))
        cursor += rendered_height + gap

    target.parent.mkdir(parents=True, exist_ok=True)
    etree.ElementTree(outer).write(
        str(target),
        encoding="utf-8",
        xml_declaration=True,
        pretty_print=False,
    )


def build_curated_supplementary_figures(
    source_dir: Path,
    target_dir: Path,
) -> dict[str, Path]:
    target_dir.mkdir(parents=True, exist_ok=True)
    image_map: dict[str, Path] = {}
    for output_number, source_numbers in SUPPLEMENTARY_FIGURE_GROUPS.items():
        for extension, composer in (
            ("png", _compose_png),
            ("pdf", _compose_pdf),
            ("svg", _compose_svg),
        ):
            source_paths = [
                source_dir / f"Supplementary_Figure_S{number}_V4.{extension}"
                for number in source_numbers
            ]
            missing = [path for path in source_paths if not path.is_file()]
            if missing:
                raise FileNotFoundError(
                    "Missing supplementary figure inputs:\n"
                    + "\n".join(str(path) for path in missing)
                )
            composer(source_paths, target_dir / f"Figure_S{output_number}.{extension}")
        image_map[f"FIGURE_S{output_number}"] = (
            target_dir / f"Figure_S{output_number}.png"
        )
    return image_map
