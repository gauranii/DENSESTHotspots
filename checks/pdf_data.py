# -----------------------------------------------------------------------------
# File:    checks/pdf_data.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Read the marks drawn in an R pdf() file (rectangles, paths, text)
#          so a figure can be compared with another by its data. Also reads
#          the PDF files matplotlib writes without compression (the figures
#          of python/zinull/plots.py), whose point markers are stamps of a
#          form drawn once.
#          Standard library only.
# Usage:   python3 checks/pdf_data.py figure.pdf   (prints a summary)
# -----------------------------------------------------------------------------
import math
import re
import sys
import zlib

TOKEN = re.compile(r'\((?:\\.|[^\\()])*\)|\[[^\]]*\]|/[^\s/\[\]()<>]+|[-+]?\d*\.?\d+|[A-Za-z\'"*]+')
PAINT = {'S': (False, True), 's': (False, True), 'f': (True, False), 'F': (True, False),
         'f*': (True, False), 'B': (True, True), 'B*': (True, True), 'b': (True, True),
         'b*': (True, True), 'n': (False, False)}


def content(path):
    """Concatenated, decompressed page content streams of an R pdf file (or
    of a matplotlib one; the forms of its point markers are left out)."""
    raw = open(path, 'rb').read()
    parts = []
    for m in re.finditer(rb'<<([^<>]*)>>\s*stream\r?\n', raw):
        if b'/Alternate' in m.group(1) or b'/FontFile' in m.group(1) or b'/Form' in m.group(1):
            continue  # colour profile, font or marker form, not page content
        body = raw[m.end():raw.index(b'endstream', m.end())]
        if b'/FlateDecode' in m.group(1):
            body = zlib.decompressobj().decompress(body)
        parts.append(body.decode('latin1'))
    return '\n'.join(parts)


def marks(path):
    """Marks of a figure in page coordinates (pt), in drawing order:
    'rects' [(x, y, w, h, fill, stroke)]: painted rectangles; fill and stroke
            are colour strings such as '0.745 0.745 0.745', None if not used.
    'paths' [{'points': [(x, y)], 'curves': n, 'closed': bool, 'fill': colour,
              'stroke': colour}]: every other painted path.
    'clips' [(x, y, w, h)]: rectangles used as clipping regions.
    'texts' [(x, y, size, string, font, rotated)].
    'stamps' [(name, x, y)]: forms drawn at a point (the point markers of a
            matplotlib figure; none in an R figure)."""
    rects, paths, clips, texts, stamps = [], [], [], [], []
    state = {'fill': '0 0 0', 'stroke': '0 0 0', 'origin': (0.0, 0.0)}
    saved, stack = [], []
    points, curves, closed, rect = [], 0, False, None
    font, tm = None, (1, 0, 0, 1, 0, 0)
    for token in TOKEN.findall(content(path)):
        if token[0] in '([/' or token[0].isdigit() or token[0] in '+-.':
            stack.append(token)
            continue
        numbers = []
        for operand in stack:
            try:
                numbers.append(float(operand))
            except ValueError:
                pass
        if token == 'q':
            saved.append(dict(state))
        elif token == 'Q' and saved:
            state = saved.pop()
        elif token in ('scn', 'sc', 'rg') and len(numbers) == 3:
            state['fill'] = ' '.join(stack[-3:])
        elif token in ('SCN', 'SC', 'RG') and len(numbers) == 3:
            state['stroke'] = ' '.join(stack[-3:])
        elif token == 'g' and numbers:
            state['fill'] = ' '.join(stack[-1:] * 3)
        elif token == 'G' and numbers:
            state['stroke'] = ' '.join(stack[-1:] * 3)
        elif token in ('m', 'l') and len(numbers) >= 2:
            points.append((numbers[-2], numbers[-1]))
        elif token == 'c' and len(numbers) >= 6:
            points.append((numbers[-2], numbers[-1]))
            curves += 1
        elif token == 'h':
            closed = True
        elif token == 're' and len(numbers) >= 4:
            rect = tuple(numbers[-4:])
        elif token in PAINT:
            fill, stroke = PAINT[token]
            fill = state['fill'] if fill else None
            stroke = state['stroke'] if stroke else None
            if rect and token == 'n':
                clips.append(rect)
            elif rect:
                rects.append(rect + (fill, stroke))
            elif points and token != 'n':
                paths.append({'points': points, 'curves': curves,
                              'closed': closed or token in ('s', 'b', 'b*'),
                              'fill': fill, 'stroke': stroke})
            points, curves, closed, rect = [], 0, False, None
        elif token == 'cm' and len(numbers) >= 6:  # a shift of the origin
            state['origin'] = (state['origin'][0] + numbers[-2], state['origin'][1] + numbers[-1])
        elif token == 'Do' and stack:
            stamps.append((stack[-1],) + state['origin'])
        elif token == 'Tf' and stack:
            font = stack[0]
        elif token == 'Tm' and len(numbers) >= 6:
            tm = tuple(numbers[-6:])
        elif token == 'Tj' and stack and stack[-1][0] == '(':
            size = math.hypot(tm[0], tm[1])
            string = re.sub(r'\\(.)', r'\1', stack[-1][1:-1])
            texts.append((tm[4], tm[5], size, string, font, abs(tm[1]) > abs(tm[0])))
        stack = []
    return {'rects': rects, 'paths': paths, 'clips': clips, 'texts': texts, 'stamps': stamps}


if __name__ == '__main__':
    for p in sys.argv[1:]:
        mk = marks(p)
        print(p, *(f'{k} {len(v)}' for k, v in mk.items()))
