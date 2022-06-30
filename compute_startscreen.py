import matplotlib.pyplot as plt

plt.axis('off')
plt.text(0.5, 0.5, "Welcome to the soliton automata software! \n Please specify your molecule below:", ha='center', va='center')
plt.savefig(f'database/startscreen.jpg', bbox_inches='tight', format='jpg', dpi=1200)
plt.show()
