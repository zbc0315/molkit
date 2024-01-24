import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, QPushButton
from PyQt5.QtGui import QPixmap
from PIL.ImageQt import ImageQt

from molkit.mol_utils import smiles_to_2d_mol, smiles_to_3d_mol


if __name__ == "__main__":
    try:
        app = QApplication(sys.argv)
        win = QMainWindow()
        win.setGeometry(400, 400, 400, 300)
        win.setWindowTitle("Smiles to 3D Mol")
        label = QLabel("Input Smiles, Please!", win)
        label.move(50, 50)
        label.resize(200, 30)
        input_smiles = QLineEdit(win)
        input_smiles.move(50, 100)
        input_smiles.resize(200, 30)

        mol_2d = QLabel(win)
        mol_2d.move(50, 150)
        mol_2d.setStyleSheet("border: 1px solid black;")
        mol_2d.resize(300, 300)

        def get_2d_mol():
            mol_2d.clear()
            smiles = input_smiles.text()
            mol_image, w, h = smiles_to_2d_mol(smiles)
            if mol_image is None:
                mol_2d.setText("Invalid Smiles!")
            else:
                print(f"w: {w}, h: {h}")
                pixmap = QPixmap.fromImage(ImageQt(mol_image))
                mol_2d.setPixmap(pixmap)
                # mol_2d.resize(int((w+1)*200), int((h+1)*200))

        input_smiles.textChanged.connect(get_2d_mol)

        button_3d = QPushButton("Convert to 3D", win)
        button_3d.move(50, 470)
        button_3d.resize(100, 30)

        def draw_3d_mol():
            smiles = input_smiles.text()
            mol_image = smiles_to_3d_mol(smiles)
            if mol_image is None:
                mol_2d.setText("Invalid Smiles!")
            else:
                print(mol_image)
                pixmap = QPixmap.fromImage(ImageQt(mol_image))
                mol_2d.setPixmap(pixmap)

        button_3d.clicked.connect(draw_3d_mol)

        win.show()
        sys.exit(app.exec_())
    except Exception as e:
        print(e)
        sys.exit(1)
