"""Tutorial images from the screenshots of the real interface (screenshots/,
made by test/screenshots_test.dart) and the scenes of
assets/tutorial/scenes.json (one per "what to do" item: screenshot, area to
point at, short label):

- assets/tutorial/shots/*.png  the screenshots the app's tutorial uses (the
  app draws the highlight and an arrow from each item to its area);
- docs/assets/tutorial/step<n>.<lang>.gif  annotated GIFs for the docs case
  study: the area highlighted, a numbered label and an arrow (straight
  lines, right angles) from the label to the area, one frame per item.

    python3 tool/make_tutorial_gifs.py      (needs Pillow)
"""
import json
import math
import pathlib
import shutil

from PIL import Image, ImageDraw, ImageFont

ROOT = pathlib.Path(__file__).resolve().parents[1]
SHOTS = ROOT / "screenshots"
APP = ROOT / "assets" / "tutorial"
DOCS = ROOT.parent / "docs" / "assets" / "tutorial"
FONT = ROOT / "assets" / "fonts" / "Montserrat-SemiBold.ttf"
WIDTH = 960
HOLD_MS = 2800
RED = (138, 0, 31)
WHITE = (255, 255, 255)

scenes = json.loads((APP / "scenes.json").read_text(encoding="utf-8"))
W0, H0 = scenes["size"]
S = WIDTH / W0


def annotate(img, rect, number, label):
    """[img] with the area dimmed around [rect], a red frame, a numbered
    label and an arrow (straight lines) from the label to the frame."""
    base = img.convert("RGBA")
    if rect is None:
        return base.convert("RGB")
    x, y, w, h = [v * S for v in rect]
    shade = Image.new("RGBA", base.size, (15, 17, 30, 95))
    ImageDraw.Draw(shade).rectangle([x, y, x + w, y + h], fill=(0, 0, 0, 0))
    out = Image.alpha_composite(base, shade)
    d = ImageDraw.Draw(out)
    d.rounded_rectangle([x - 3, y - 3, x + w + 3, y + h + 3], radius=8, outline=RED, width=4)

    font = ImageFont.truetype(str(FONT), 17)
    tw = d.textlength(label, font=font)
    lw, lh = tw + 58, 40
    W, H = out.size
    cx, cy = x + w / 2, y + h / 2
    # The label on the side away from the area, inside the picture.
    lx = W * 0.56 if cx < W / 2 else W * 0.05
    ly = H * 0.72 if cy < H / 2 else H * 0.14
    lx = min(max(lx, 10), W - lw - 10)
    if lx < x + w and lx + lw > x and ly < y + h and ly + lh > y:
        ly = max(10, y - lh - 40) if y > H / 2 else min(H - lh - 10, y + h + 40)
    # Arrow in straight lines with a right angle: from the side of the
    # label facing the area, sideways to the area's x, then up or down into
    # it; the label is drawn over the line.
    ex = min(max(lx + lw / 2, x + 12), x + w - 12)
    below = ly > y + h
    above = ly + lh < y
    if below or above:
        end = (ex, y + h + 4 if below else y - 4)
        if lx <= ex <= lx + lw:
            start = (ex, ly if below else ly + lh)
            pts = [start, end]
        else:
            start = (lx if ex < lx else lx + lw, ly + lh / 2)
            pts = [start, (ex, start[1]), end]
    else:
        # Beside the area: straight across to its nearest side.
        sy = min(max(ly + lh / 2, y + 6), y + h - 6)
        start = (lx if x + w < lx else lx + lw, sy)
        end = (x + w + 4 if x + w < lx else x - 4, sy)
        pts = [start, end]
    d.line(pts, fill=RED, width=4, joint="curve")
    (ax, ay), (bx, by) = pts[-2], pts[-1]
    a = math.atan2(by - ay, bx - ax)
    head = [(bx, by), (bx - 16 * math.cos(a - 0.45), by - 16 * math.sin(a - 0.45)), (bx - 16 * math.cos(a + 0.45), by - 16 * math.sin(a + 0.45))]
    d.polygon(head, fill=RED)
    d.rounded_rectangle([lx, ly, lx + lw, ly + lh], radius=10, fill=WHITE, outline=RED, width=3)
    d.ellipse([lx + 8, ly + 8, lx + 32, ly + 32], fill=RED)
    num = str(number)
    d.text((lx + 20 - d.textlength(num, font=font) / 2, ly + 9), num, font=font, fill=WHITE)
    d.text((lx + 42, ly + 9), label, font=font, fill=(31, 41, 55))
    return out.convert("RGB")


def load(name):
    im = Image.open(SHOTS / f"{name}.png").convert("RGB")
    return im.resize((WIDTH, round(im.height * S)), Image.LANCZOS)


# The app: the screenshots as they are.
(APP / "shots").mkdir(parents=True, exist_ok=True)
for name in sorted({sc["shot"] for page in scenes["pages"] for sc in page}):
    shutil.copy(SHOTS / f"{name}.png", APP / "shots" / f"{name}.png")

# The docs: one annotated GIF per step and language.
DOCS.mkdir(parents=True, exist_ok=True)
for old in DOCS.glob("step*.gif"):
    old.unlink()
for step, page in enumerate(scenes["pages"]):
    for lang in ("en", "pt"):
        frames = [annotate(load(sc["shot"]), sc["rect"], i + 1, sc[lang]) for i, sc in enumerate(page)]
        pal = [f.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE) for f in frames]
        target = DOCS / f"step{step}.{lang}.gif"
        pal[0].save(target, save_all=True, append_images=pal[1:], duration=HOLD_MS, loop=0, optimize=True, disposal=1)
        print(target.relative_to(ROOT.parent), f"{target.stat().st_size // 1024} KB", len(frames), "frames")
