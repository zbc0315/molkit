import tkinter as tk
from PIL import ImageTk, Image
from molkit.mol_utils import smiles_to_2d_mol, smiles_to_3d_mol


if __name__ == "__main__":
    font = ('Arial', 20)
    root = tk.Tk()
    label_smiles = tk.Label(root, text="Input Smiles, Please!", font=font)
    label_smiles.pack()
    entry_smiles = tk.Entry(root, width=20, font=font)
    entry_smiles.pack(ipady=10)

    label_mol = tk.Label(root)
    def draw_mol(event):
        label_mol.config(image=None)
        smiles = event.widget.get()
        mol_image, w, h = smiles_to_2d_mol(smiles)
        if mol_image is not None:
            tk_img = ImageTk.PhotoImage(mol_image)
            label_mol.config(image=tk_img)
            label_mol.pack()
        else:
            label_mol.config(text="Invalid Smiles!")
            label_mol.pack()

    entry_smiles.bind("<Key>", draw_mol)

    root.mainloop()

    # root = Tk()
    # frm = ttk.Frame(root, padding=10)
    # frm.grid()
    # ttk.Label(frm, text="Hello World!").grid(column=0, row=0)
    # ttk.Entry(frm, width=20).grid(column=0, row=1)
    # ttk.Button(frm, text="Quit", command=root.destroy).grid(column=0, row=2)
    # root.mainloop()
