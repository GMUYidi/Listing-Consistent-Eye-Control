"""Add a 4 pt safety margin to Figure 5 without changing its vector artwork."""
from pathlib import Path
import shutil
import sys

from pypdf import PdfReader, PdfWriter
from pypdf.generic import RectangleObject

def main():
    pdf = (Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else
           Path(__file__).resolve().parents[2] / 'outputs/figure5/net_torque_eye_orientation.pdf')
    reader = PdfReader(pdf)
    if (reader.metadata or {}).get('/Figure5MarginPoints') == '4':
        print('Figure 5 already has the 4 pt safety margin.')
        return
    raw_pdf = pdf.with_name(pdf.stem + '.matlab_raw.pdf')
    shutil.copyfile(pdf, raw_pdf)
    writer = PdfWriter()
    writer.clone_document_from_reader(reader)
    for page in writer.pages:
        box = page.mediabox
        padded = RectangleObject([float(box.left)-4, float(box.bottom)-4,
                                   float(box.right)+4, float(box.top)+4])
        page.mediabox = padded
        page.cropbox = padded
    writer.add_metadata({'/Figure5MarginPoints':'4'})
    with pdf.open('wb') as stream:
        writer.write(stream)
    print('Figure 5 PDF margin added; original vector artwork preserved.')

if __name__ == '__main__':
    main()
