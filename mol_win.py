import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import ImageTk, Image


from rdkit.Chem import AllChem, Draw


class MolWin(tk.Tk):

    def __init__(self):
        super().__init__()
        self.mol_3d = None
        self.title("MolWin")
        self.geometry("800x600")
        self.resizable(False, False)
        self.font = ('Arial', 20)
        self.label_smiles = tk.Label(self, text="Input Smiles, Please!", font=self.font)
        self.label_smiles.pack()
        self.entry_smiles = tk.Entry(self, width=20, font=self.font)
        self.entry_smiles.pack(ipady=10)

        self.label_mol = tk.Label(self)
        self.tk_image = None

        self.entry_smiles.bind("<Return>", self.draw_mol)

        self.button_to_3d = tk.Button(self, text="To 3D", font=self.font, command=self.draw_3d_mol)
        self.button_to_3d.pack()

        self.button_save = tk.Button(self, text="Save", font=self.font, command=self.save_file)
        self.button_save.pack()

    def save_file(self):
        if self.mol_3d is None:
            messagebox.showinfo("Alert", "'To 3D' first, Please!")
        else:
            file_path = filedialog.asksaveasfilename(defaultextension=".mol",
                                                     filetypes=[("Text Files", "*.mol"), ("All Files", "*.*")])

            if file_path:
                # 在此处编写保存文件的逻辑
                self.write_mol(file_path)
            else:
                messagebox.showinfo("Alert", "Please input file name!")

    def draw_mol(self, event):
        smiles = event.widget.get()
        mol_image, w, h = self.smiles_to_2d_mol(smiles)
        if mol_image is not None:
            tk_img = ImageTk.PhotoImage(mol_image)
            self.label_mol.config(image=tk_img)
            self.label_mol.image = tk_img
            self.label_mol.pack()
        else:
            print("Invalid Smiles!")
            self.label_mol.config(image=None)
            self.label_mol.image = None
            self.label_mol.pack()

    def draw_3d_mol(self):
        # nonlocal tk_image, smiles
        smiles = self.entry_smiles.get()
        mol_image = self.smiles_to_3d_mol(smiles)
        print(smiles)
        print(f"3D mol_image: {mol_image}")
        if mol_image is not None:
            tk_img = ImageTk.PhotoImage(mol_image)
            self.label_mol.config(image=tk_img)
            self.label_mol.image = tk_img
            self.label_mol.pack()
        else:
            self.label_mol.config(image=None)
            self.label_mol.image = None
            # self.label_mol.config(text="Invalid Smiles!")
            # self.label_mol.text = "Invalid Smiles!"
            self.label_mol.pack()

    def smiles_to_3d_mol(self, smiles):
        try:
            mol = AllChem.MolFromSmiles(smiles)
            mol = AllChem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol)
            self.mol_3d = mol
            return Draw.MolToImage(mol)
        except:
            return None

    def write_mol(self, fp):
        with open(fp, "w", encoding="utf-8") as f:
            f.write(AllChem.MolToMolBlock(self.mol_3d))
        messagebox.showinfo("Alert", "Save Success to {}".format(fp))

    def smiles_to_2d_mol(self, smiles):
        try:
            mol = AllChem.MolFromSmiles(smiles)
            AllChem.Compute2DCoords(mol)
            coords = mol.GetConformer().GetPositions()
            max_x = max(coords, key=lambda x: x[0])[0]
            min_x = min(coords, key=lambda x: x[0])[0]
            max_y = max(coords, key=lambda x: x[1])[1]
            min_y = min(coords, key=lambda x: x[1])[1]
            return Draw.MolToImage(mol), max_x - min_x, max_y - min_y
        except:
            return None, None, None


if __name__ == "__main__":
    mw = MolWin()
    mw.mainloop()
