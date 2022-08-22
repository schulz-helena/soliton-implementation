import io

import matplotlib.pyplot as plt
from PIL import Image


class Startscreen:
    """Computes image for the start screen.
    """

    def __init__(self):
        plt.axis('off')
        plt.text(0.5, 0.5, "Welcome to the Mini Soliton Automata Software! \n Please specify your molecule below:", ha='center', va='center', fontname = "Futura", fontsize = 12)
        buf = io.BytesIO()
        plt.savefig(buf, bbox_inches='tight', format='jpg', dpi=800)
        buf.seek(0)
        self.image = Image.open(buf)
        self.image: Image = self.image.convert("RGBA")
        """Welcoming image."""
        buf.close()
