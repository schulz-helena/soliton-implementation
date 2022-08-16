import io

import matplotlib.pyplot as plt
from PIL import Image


class Startscreen:

    def __init__(self):
        plt.axis('off')
        plt.text(0.5, 0.5, "Welcome to the soliton automata software! \n Please specify your molecule below:", ha='center', va='center', fontname = "Futura", fontsize = 12)
        #plt.savefig(f'database/startscreen.jpg', bbox_inches='tight', format='jpg', dpi=1200)
        #plt.show()
        buf = io.BytesIO()
        plt.savefig(buf, bbox_inches='tight', format='jpg', dpi=800)
        buf.seek(0)
        self.image = Image.open(buf)
        self.image = self.image.convert("RGBA")
        buf.close()
