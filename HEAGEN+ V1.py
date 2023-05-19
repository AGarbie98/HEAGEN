import random
import tkinter as tk
from tkinter import simpledialog, messagebox
from tkinter import messagebox
from ase import Atoms
from ase.spacegroup import crystal
import numpy as np


ELEMENTS = [
    "Li", "Be", "Na", "Mg", "K", "Ca", "Rb", "Sr", "Cs", "Ba", "La", "Hf",
    "Y", "Sc", "Ti", "Zr", "V", "Nb", "Ta", "Cr", "Mo", "W", "Mn", "Tc",
    "Re", "Fe", "Ru", "Os", "Co", "Rh", "Ir", "Ni", "Pd", "Pt", "Cu", "Ag",
    "Au", "Zn", "Cd", "Hg", "B", "Al", "Ga", "In", "Tl", "C", "Si", "Ge",
    "Sn", "Pb", "N", "P", "As", "Sb", "Bi"
]

LATTICE_CONSTANTS = {
    "Li": 3.49, "Be": 2.29, "Na": 4.23, "Mg": 3.21, "K": 5.23, "Ca": 5.58,
    "Sc": 3.31, "Ti": 2.95, "V": 3.03, "Cr": 2.88, "Mn": 8.89, "Fe": 2.87,
    "Co": 2.51, "Ni": 3.52, "Cu": 3.61, "Zn": 2.66,"Ga": 4.51, "Ge": 5.66,
    "As": 4.13, "Rb": 5.59, "Sr": 6.08, "Y": 3.65, "Zr": 3.23,"Nb": 3.30,
    "Mo": 3.15, "Tc": 2.74, "Ru": 2.70, "Rh": 3.80, "Pd": 3.89, "Ag": 4.09,
    "Cd": 2.98, "In": 4.59,"Sn": 5.83, "Sb": 4.31, "Cs": 6.05, "Ba": 5.02,
    "La": 3.75, "Hf": 3.20, "Ta": 3.31, "W": 3.16, "Re": 2.76, "Os": 2.74,
    "Ir": 3.84,"Pt": 3.92,"Au": 4.08,"Hg": 3.46,"Tl": 3.46, "Pb": 4.95,
    "Bi": 4.75,  "Al": 4.05,  "Si": 5.43,  "B": 8.73,  "C": 3.57,  "N": 4.039,
    "P": 7.17,
}

class ElementSelector(tk.Tk):
    def __init__(self, elements):
        super().__init__()
        self.title("Element Selector")
        self.geometry("500x700")
        self.elements = elements
        self.selected_elements = []
        self.buttons = {}
        self.periodic_table_layout = [
            ['Li', 'Be', '', '', '', '', '',' ', '', '', '', '', 'B', 'C', 'N'],
            ['Na', 'Mg', '', '', '', '', '', '', '', '', '', '', 'Al', 'Si', 'P'],
            ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As'],
            ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb'],
            ['Cs', 'Ba', 'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi'],
        ]

        for i, row in enumerate(self.periodic_table_layout):
            for j, element in enumerate(row):
                if element in elements:
                    button = tk.Button(self, text=element, width=5, height=3, command=lambda e=element: self.select_element(e))
                    button.grid(row=i, column=j)
                    self.buttons[element] = button

        self.submit_button = tk.Button(self, text="Submit", command=self.submit)
        self.submit_button.grid(row=len(self.periodic_table_layout) + 1, column=0, columnspan=len(self.periodic_table_layout[0]))

    def select_element(self, element):
        if element in self.selected_elements:
            self.selected_elements.remove(element)
            self.buttons[element].config(relief=tk.RAISED)
        else:
            self.selected_elements.append(element)
            self.buttons[element].config(relief=tk.SUNKEN)

    def submit(self):
        if len(self.selected_elements) < 1:
            messagebox.showwarning("Warning", "Please select at least one element.")
        else:
            self.quit()

    def run_and_get_elements(self):
        self.mainloop()
        return self.selected_elements

###############################

class CrystalStructureSelector(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Crystal Structure Selector")
        self.geometry("200x100")
        self.crystal_structure = None

        self.fcc_button = tk.Button(self, text="FCC", command=lambda: self.select_structure("FCC"))
        self.fcc_button.pack(side=tk.LEFT, padx=20)

        self.bcc_button = tk.Button(self, text="BCC", command=lambda: self.select_structure("BCC"))
        self.bcc_button.pack(side=tk.LEFT, padx=20)

    def select_structure(self, structure):
        self.crystal_structure = structure
        self.quit()

    def run_and_get_structure(self):
        self.mainloop()
        return self.crystal_structure

def get_elements():
    selector = ElementSelector(ELEMENTS)
    selected_elements = selector.run_and_get_elements()
    return selected_elements

def get_crystal_structure():
    selector = CrystalStructureSelector()
    crystal_structure = selector.run_and_get_structure()
    return crystal_structure

def get_elements():
    selector = ElementSelector(ELEMENTS)
    selected_elements = selector.run_and_get_elements()
    return selected_elements

########################

class MolarRatioSelector(tk.Tk):
    def __init__(self, elements):
        super().__init__()
        self.title("Molar Ratio Selector")
        self.geometry("300x" + str(30 * len(elements) + 50))
        self.elements = elements
        self.ratios = {}
        self.entries = {}
        for i, element in enumerate(elements):
            tk.Label(self, text=element).grid(row=i, column=0)
            entry = tk.Entry(self)
            entry.grid(row=i, column=1)
            self.entries[element] = entry

        self.submit_button = tk.Button(self, text="Submit", command=self.submit)
        self.submit_button.grid(row=len(self.elements), column=0, columnspan=2)

    def submit(self):
        for element, entry in self.entries.items():
            try:
                self.ratios[element] = int(entry.get())
            except ValueError:
                messagebox.showwarning("Warning", "Please enter a valid integer for each element.")
                return
        self.quit()

    def run_and_get_ratios(self):
        self.mainloop()
        return self.ratios

##########################
def calculate_HEA_lattice_constant(elements, molar_ratios):
    total_ratio = sum(molar_ratios.values())
    lattice_constant = 0
    for element in elements:
        ratio = molar_ratios[element]
        element_lattice_constant = LATTICE_CONSTANTS[element]
        lattice_constant += ratio / total_ratio * element_lattice_constant
    return lattice_constant

def get_molar_ratios(elements):
    selector = MolarRatioSelector(elements)
    molar_ratios = selector.run_and_get_ratios()
    return molar_ratios

def generate_HEA_supercell(elements, lattice_constant, supercell_size, crystal_structure, molar_ratios):
    # Ensure the supercell size is a multiple of the number of elements
    total_atoms = supercell_size[0] * supercell_size[1] * supercell_size[2]
    total_atoms *= 2 if crystal_structure == "BCC" else 4  # BCC has 2 atoms per unit cell, FCC has 4
    if total_atoms % sum(molar_ratios.values()) != 0:
        raise ValueError("The total number of atoms in the supercell must be a multiple of the sum of the molar ratios.")

    # Create unit cell
    basis = [(0, 0, 0), (0.5, 0.5, 0.5)] if crystal_structure == "BCC" else [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    spacegroup = 229 if crystal_structure == "BCC" else 225
    unit_cell = crystal(
        symbols=elements,
        basis=basis,
        spacegroup=spacegroup,
        cellpar=[lattice_constant, lattice_constant, lattice_constant, 90, 90, 90]
    )

    # Create supercell
    supercell = unit_cell.repeat(supercell_size)

    # Randomly distribute elements in the supercell
    symbols = []
    for element, ratio in molar_ratios.items():
        symbols += [element] * ratio
    symbols *= total_atoms // len(symbols)
    random.shuffle(symbols)
    supercell.set_chemical_symbols(symbols)

    return supercell

def main():
    # Create a Tkinter root window (it won't be displayed)
    root = tk.Tk()
    root.withdraw()

    # Get the selected elements
    elements = get_elements()

    # Get the molar ratios
    molar_ratios = get_molar_ratios(elements)

    # Calculate the lattice constant of the HEA
    lattice_constant = calculate_HEA_lattice_constant(elements, molar_ratios)

     # Display the selected elements, molar ratios, and calculated lattice constant
    info = "Selected elements: " + ", ".join(elements) + "\n"
    info += "Molar ratios: " + ", ".join(f"{element}: {ratio}" for element, ratio in molar_ratios.items()) + "\n"
    info += "Calculated lattice constant: " + str(lattice_constant)
    messagebox.showinfo("HEA Information", info)

    # Get the crystal structure
    crystal_structure = get_crystal_structure()
    
    # Calculate the smallest supercell size that is a multiple of the sum of the molar ratios
    sum_ratios = sum(molar_ratios.values())
    supercell_size = (sum_ratios, sum_ratios, sum_ratios)

    try:
        # Generate the high-entropy alloy structure
        structure = generate_HEA_supercell(elements, lattice_constant, supercell_size, crystal_structure, molar_ratios)

        # Save the structure to a .cif file
        filename = "".join(elements) + "_" + crystal_structure + ".cif"
        structure.write(filename)

        messagebox.showinfo("Success", "The structure was successfully generated and saved to " + filename)
    except Exception as e:
        messagebox.showerror("Error", "An error occurred while generating the structure: " + str(e))

if __name__ == "__main__":
    main()
