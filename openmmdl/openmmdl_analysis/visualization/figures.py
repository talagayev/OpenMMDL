import os
import re

import cairosvg
import MDAnalysis as mda
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, rdCoordGen, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import pylab
from PIL import Image


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "HID": "H", "HIE": "H", "HIP": "H",
}

INTERACTION_COLORS = {
    "hbond": "#4C78A8",
    "hydrophobic": "#F58518",
    "waterbridge": "#72B7B2",
    "saltbridge": "#E45756",
    "pistacking": "#B279A2",
    "pication": "#9D755D",
    "halogen": "#54A24B",
    "metal": "#BAB0AC",
}


class FigureMerger:
    """
    Handles the creation and merging of binding mode figures with corresponding legends.

    Attributes
    ----------
    binding_mode : str
        Name of the binding mode for which figures are created.
    occurrence_percent : float
        Occurrence percentage of the binding mode.
    split_data : list of str
        Interaction descriptors used for generating the figure legend.
    merged_image_paths : list of str
        List storing paths to the output merged images.
    """

    def __init__(self, binding_mode, occurrence_percent, split_data, merged_image_paths):
        self.binding_mode = binding_mode
        self.occurrence_percent = occurrence_percent
        self.split_data = split_data
        self.merged_image_paths = merged_image_paths

    def create_and_merge_images(self):
        """
        Create and merge images to generate a legend for binding modes.

        Returns
        -------
        list of str
            Updated list of paths to the merged images.
        """
        # Create the main figure and axis
        fig = pylab.figure()
        ax = fig.add_subplot(111)

        # Data for the x-axis and random data for demonstration
        x = range(10)
        data_points = [pylab.randn(10) for _ in range(len(self.split_data))]

        # Plot lines on the same axis and collect them into a list
        lines = []
        filtered_split_data = [entry for entry in self.split_data if "FRAME" not in entry]
        for i, data in enumerate(filtered_split_data):
            y = data_points[i]
            label = data.split()[-1]
            type = data.split()[-2]
            if label == "hydrophobic":
                (line,) = ax.plot(x, y, label=data, color=(1.0, 1.0, 0.0), linewidth=5.0)  # yellow
            elif label == "hbond":
                if type == "Acceptor":
                    (line,) = ax.plot(x, y, label=data, color=(1.0, 0.6, 0.6), linewidth=5.0)  # light red / pink
                elif type == "Donor":
                    (line,) = ax.plot(x, y, label=data, color=(0.3, 0.5, 1.0), linewidth=5.0)  # light blue
            elif label == "halogen":
                (line,) = ax.plot(x, y, label=data, color=(1.0, 0.0, 0.9), linewidth=5.0)  # magenta / hot pink
            elif label == "pistacking":
                (line,) = ax.plot(x, y, label=data, color=(0.0, 0.0, 1.0), linewidth=5.0)  # blue
            elif label == "pication":
                (line,) = ax.plot(x, y, label=data, color=(0.0, 0.0, 1.0), linewidth=5.0)  # blue
            elif label == "waterbridge":
                (line,) = ax.plot(x, y, label=data, color=(0.0, 1.0, 0.9), linewidth=5.0)  # cyan / aqua
            elif label == "metal":
                (line,) = ax.plot(x, y, label=data, color=(1.0, 0.6, 0.0), linewidth=5.0)  # orange
            elif label == "saltbridge":
                if type == "NI":
                    (line,) = ax.plot(x, y, label=data, color=(1.0, 0.6, 0.0), linewidth=5.0)  # orange
                elif type == "PI":
                    (line,) = ax.plot(x, y, label=data, color=(0.3, 0.9, 0.8), linewidth=5.0)  # turquoise / teal
            else:
                (line,) = ax.plot(x, y, label=data)
            lines.append(line)

        # Create a separate figure for the legend
        figlegend = pylab.figure(figsize=(8, 6))

        # Add a legend to the subplot (ax) using the lines and full entries as labels
        legend = figlegend.legend(lines, filtered_split_data, loc="center")

        # Set the linewidth of the legend lines to be thicker
        for line in legend.get_lines():
            line.set_linewidth(5.0)

        # Add text above the legend
        figlegend.text(0.5, 0.9, f"{self.binding_mode}", ha="center", fontsize=12, weight="bold")
        figlegend.text(
            0.5,
            0.85,
            f"Occurrence {self.occurrence_percent}%",
            ha="center",
            fontsize=12,
            weight="bold",
        )

        # Save the legend figure to a file
        legend_filename = f"{self.binding_mode}_legend.png"
        figlegend.savefig(legend_filename)

        # Read the two images
        image1 = Image.open(f"{self.binding_mode}.png")
        image2 = Image.open(legend_filename)

        # Resize the first image
        image1_size = image1.size
        image2_size = image2.size
        total_width = image1_size[0] + image2_size[0]
        new_image = Image.new("RGB", (total_width, image1_size[1]))
        new_image.paste(image1, (0, 0))
        new_image.paste(image2, (image1_size[0], 0))

        # Save the merged image
        merged_image_filename = f"{self.binding_mode}_merged.png"
        new_image.save(merged_image_filename, "PNG")

        # Append the merged image path to the list
        self.merged_image_paths.append(merged_image_filename)

        # Remove the original files
        os.remove(f"{self.binding_mode}.png")
        os.remove(legend_filename)
        os.remove(f"{self.binding_mode}.svg")

        return self.merged_image_paths


class FigureArranger:
    """
    Arranges multiple merged binding mode figures into a single image.

    Attributes
    ----------
    merged_image_paths : list of str
        List of file paths to pre merged figures.
    output_path : str
        File path where the final arranged figure should be saved.
    """

    def __init__(self, merged_image_paths, output_path):
        self.merged_image_paths = merged_image_paths
        self.output_path = output_path

    def arranged_figure_generation(self):
        """
        Generate an arranged figure by arranging merged images in rows and columns.

        Returns
        -------
        None
            This function writes out a figure and does not return anything.
        """
        # Open the list of images
        merged_images = [Image.open(path) for path in self.merged_image_paths]

        # Calculate the maximum width and height for the images
        max_width = max(image.size[0] for image in merged_images)
        max_height = max(image.size[1] for image in merged_images)

        # Determine the number of images per row (in your case, 2 images per row)
        images_per_row = 2

        # Calculate the number of rows and columns required
        num_rows = (len(merged_images) + images_per_row - 1) // images_per_row
        total_width = max_width * images_per_row
        total_height = max_height * num_rows

        # Create a new image with the calculated width and height
        big_figure = Image.new("RGB", (total_width, total_height), (255, 255, 255))  # Set background to white

        x_offset = 0
        y_offset = 0

        for image in merged_images:
            # Paste the image onto the big figure
            big_figure.paste(image, (x_offset, y_offset))

            # Update offsets
            x_offset += max_width

            # Move to the next row if necessary
            if x_offset >= total_width:
                x_offset = 0
                y_offset += max_height

        # Save the big figure
        big_figure.save(self.output_path, "PNG")

        # Ensure target directories exist
        target_dir = "Binding_Modes_Markov_States"
        individual_dir = os.path.join(target_dir, "individual_figures")
        os.makedirs(target_dir, exist_ok=True)
        os.makedirs(individual_dir, exist_ok=True)

        # Move the arranged overview figure into the main binding-mode directory
        new_path = os.path.join(target_dir, os.path.basename(self.output_path))
        os.replace(self.output_path, new_path)

        # Move the individual binding-mode figures into a dedicated subfolder
        for path in self.merged_image_paths:
            individual_path = os.path.join(individual_dir, os.path.basename(path))
            if os.path.abspath(path) != os.path.abspath(individual_path):
                os.replace(path, individual_path)


class PeptideBindingModeFigureGenerator:
    def __init__(self, complex_pdb_file, peptide_chain_id, rdkit_max_residues=12):
        self.complex_pdb_file = complex_pdb_file
        self.peptide_chain_id = peptide_chain_id
        self.rdkit_max_residues = rdkit_max_residues

    def _load_peptide_residues(self):
        u = mda.Universe(self.complex_pdb_file)
        peptide_atoms = u.select_atoms(f"chainID {self.peptide_chain_id}")

        residues = []
        for residue in peptide_atoms.residues:
            resname3 = str(residue.resname).upper()
            residues.append(
                {
                    "resnum": int(residue.resid),
                    "resname3": resname3,
                    "one_letter": AA3_TO_1.get(resname3, "X"),
                }
            )

        return residues

    def _extract_interaction_annotations(self, interaction_values):
        """
        Parse peptide-mode binding-mode columns like:
          44GLUA_176LYS_LYS_PI_saltbridge
          81ASPA_182LYS_LYS_PI_saltbridge
          53SERA_180TYR_hydrophobic
          12GLUA_177LYS_Donor_hbond

        Returns
        -------
        residue_to_interactions : dict[int, set[str]]
            residue number -> interaction types
        residue_to_details : dict[int, list[dict]]
            residue number -> detailed interaction descriptions
        """
        residue_to_interactions = {}
        residue_to_details = {}

        for value in interaction_values:
            parts = str(value).split("_")
            if len(parts) < 3:
                continue

            protein_token = parts[0]
            peptide_token = parts[1]
            interaction_type = parts[-1].lower()
            extra_detail = "_".join(parts[2:-1]) if len(parts) > 3 else ""

            match = re.match(r"(\d+)([A-Za-z0-9]+)", peptide_token)
            if match is None:
                continue

            resnum = int(match.group(1))
            peptide_label = peptide_token

            residue_to_interactions.setdefault(resnum, set()).add(interaction_type)
            residue_to_details.setdefault(resnum, [])

            entry = {
                "peptide_label": peptide_label,
                "protein_label": protein_token,
                "interaction_type": interaction_type,
                "extra_detail": extra_detail,
            }

            if entry not in residue_to_details[resnum]:
                residue_to_details[resnum].append(entry)

        return residue_to_interactions, residue_to_details

    def _hex_to_rgb_tuple(self, hex_color):
        hex_color = hex_color.lstrip("#")
        return tuple(int(hex_color[i:i + 2], 16) / 255 for i in (0, 2, 4))

    def _primary_interaction_type(self, interaction_types):
        priority = [
            "saltbridge",
            "hbond",
            "pication",
            "pistacking",
            "hydrophobic",
            "waterbridge",
            "halogen",
            "metal",
        ]

        for interaction_type in priority:
            if interaction_type in interaction_types:
                return interaction_type

        return sorted(interaction_types)[0]

    def _interaction_color_rgb(self, interaction_type):
        return self._hex_to_rgb_tuple(
            INTERACTION_COLORS.get(interaction_type, "#444444")
        )

    def _write_peptide_only_pdb(self, output_pdb):
        u = mda.Universe(self.complex_pdb_file)
        peptide_atoms = u.select_atoms(f"chainID {self.peptide_chain_id}")
        peptide_atoms.write(output_pdb)
        return output_pdb

    def _draw_rdkit_peptide(
        self,
        peptide_pdb,
        residue_to_interactions,
        residue_label_map,
        output_svg,
    ):
        mol = Chem.MolFromPDBFile(
            peptide_pdb,
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
        if mol is None:
            raise ValueError(f"Could not create RDKit molecule from {peptide_pdb}")

        Chem.SanitizeMol(mol, catchErrors=True)
        mol = Chem.RemoveHs(mol)

        # Prefer CoordGen for cleaner 2D peptide layouts; fall back to rdDepictor.
        try:
            rdCoordGen.AddCoords(mol)
        except Exception:
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(mol, canonOrient=True, clearConfs=True)

        highlight_atoms = []
        highlight_atom_colors = {}
        note_atom_by_residue = {}

        for atom in mol.GetAtoms():
            residue_info = atom.GetPDBResidueInfo()
            if residue_info is None:
                continue

            resnum = residue_info.GetResidueNumber()
            if resnum in residue_to_interactions:
                atom_idx = atom.GetIdx()
                highlight_atoms.append(atom_idx)
                primary_interaction = self._primary_interaction_type(
                    residue_to_interactions[resnum]
                )
                highlight_atom_colors[atom_idx] = self._interaction_color_rgb(
                    primary_interaction
                )

                atom_name = residue_info.GetName().strip()
                if resnum not in note_atom_by_residue or atom_name == "CA":
                    note_atom_by_residue[resnum] = atom_idx

        for resnum, label_id in residue_label_map.items():
            atom_idx = note_atom_by_residue.get(resnum)
            if atom_idx is not None:
                mol.GetAtomWithIdx(atom_idx).SetProp("atomNote", str(label_id))

        drawer = rdMolDraw2D.MolDraw2DSVG(1400, 500)
        options = drawer.drawOptions()
        options.padding = 0.08
        options.fixedBondLength = 32
        options.additionalAtomLabelPadding = 0.15
        options.minFontSize = 10
        options.maxFontSize = 18
        options.annotationFontScale = 1.2

        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=highlight_atom_colors,
        )
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace("svg:", "")
        with open(output_svg, "w") as f:
            f.write(svg)

    def generate_binding_mode_figure(
        self,
        binding_mode,
        interaction_values,
        occurrence_percent,
        output_dir=".",
    ):
        residues = self._load_peptide_residues()
        residue_to_interactions, residue_to_details = self._extract_interaction_annotations(
            interaction_values
        )
        residue_label_map = self._build_residue_label_map(residue_to_details)
        title = f"{binding_mode} ({occurrence_percent:.1f}%)"

        if len(residues) <= self.rdkit_max_residues:
            peptide_pdb = os.path.join(output_dir, f"{binding_mode}_peptide_only.pdb")
            svg_path = os.path.join(output_dir, f"{binding_mode}.svg")
            main_png_path = os.path.join(output_dir, f"{binding_mode}_main.png")
            legend_png_path = os.path.join(output_dir, f"{binding_mode}_legend.png")
            merged_png_path = os.path.join(output_dir, f"{binding_mode}.png")

            self._write_peptide_only_pdb(peptide_pdb)
            self._draw_rdkit_peptide(
                peptide_pdb,
                residue_to_interactions,
                residue_label_map,
                svg_path,
            )
            cairosvg.svg2png(url=svg_path, write_to=main_png_path)
            self._draw_interaction_legend(
                residue_label_map,
                residue_to_details,
                legend_png_path,
                title=title,
            )
            self._merge_main_and_legend(
                main_png_path,
                legend_png_path,
                merged_png_path,
            )

            return merged_png_path

        main_png_path = os.path.join(output_dir, f"{binding_mode}_main.png")
        legend_png_path = os.path.join(output_dir, f"{binding_mode}_legend.png")
        merged_png_path = os.path.join(output_dir, f"{binding_mode}.png")

        self._draw_sequence_diagram(
            residues,
            residue_to_interactions,
            residue_label_map,
            main_png_path,
            title=title,
        )

        self._draw_interaction_legend(
            residue_label_map,
            residue_to_details,
            legend_png_path,
            title=title,
        )
        self._merge_main_and_legend(
            main_png_path,
            legend_png_path,
            merged_png_path,
        )

        return merged_png_path

    def _draw_sequence_diagram(
        self,
        residues,
        residue_to_interactions,
        residue_label_map,
        output_png,
        title=None,
    ):
        n_residues = len(residues)
        fig_width = max(12, n_residues * 0.6)

        fig, ax = plt.subplots(figsize=(fig_width, 3.8))
        ax.hlines(y=0, xmin=0, xmax=max(n_residues - 1, 0), linewidth=1.5)

        for i, residue in enumerate(residues):
            ax.text(
                i,
                0,
                residue["one_letter"],
                ha="center",
                va="center",
                fontsize=12,
            )
            ax.text(
                i,
                -0.45,
                str(residue["resnum"]),
                ha="center",
                va="center",
                fontsize=8,
            )

            labels = residue_to_interactions.get(residue["resnum"], set())
            if not labels:
                continue

            interaction_type = self._primary_interaction_type(labels)
            color = INTERACTION_COLORS.get(interaction_type, "#444444")

            ax.annotate(
                "",
                xy=(i, 0.18),
                xytext=(i, 0.95),
                arrowprops=dict(
                    arrowstyle="->",
                    lw=1.2,
                    color=color,
                ),
            )
            ax.text(
                i,
                1.15,
                str(residue_label_map[residue["resnum"]]),
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                color=color,
            )

        ax.text(-0.8, 0, "N-term", ha="right", va="center", fontsize=10)
        ax.text(n_residues - 1 + 0.8, 0, "C-term", ha="left", va="center", fontsize=10)

        if title:
            ax.set_title(title)

        ax.set_xlim(-1.5, n_residues + 0.5)
        ax.set_ylim(-0.9, 1.9)
        ax.axis("off")

        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _build_residue_label_map(self, residue_to_details):
        """
        Assign compact integer labels (1, 2, 3, ...) to interacting peptide residues.
        """
        interacting_residues = sorted(residue_to_details.keys())
        return {resnum: i + 1 for i, resnum in enumerate(interacting_residues)}

    def _draw_interaction_legend(self, residue_label_map, residue_to_details, output_png, title=None):
        """
        Draw a side legend like:
          1. 176LYS
             - saltbridge ↔ 44GLUA (LYS_PI)
             - hbond ↔ 12GLUA (Donor)
        """
        line_count = 2
        for resnum in sorted(residue_to_details):
            line_count += 1 + len(residue_to_details[resnum])

        fig_height = max(4, line_count * 0.35)
        fig, ax = plt.subplots(figsize=(7.5, fig_height))
        ax.axis("off")

        y = 0.98
        if title:
            ax.text(0.00, y, title, fontsize=12, fontweight="bold", va="top", transform=ax.transAxes)
            y -= 0.08

        for resnum in sorted(residue_to_details):
            label_id = residue_label_map[resnum]
            entries = residue_to_details[resnum]
            peptide_label = entries[0]["peptide_label"]

            ax.text(
                0.00,
                y,
                f"{label_id}. {peptide_label}",
                fontsize=11,
                fontweight="bold",
                va="top",
                transform=ax.transAxes,
            )
            y -= 0.06

            for entry in entries:
                interaction_type = entry["interaction_type"]
                protein_label = entry["protein_label"]
                extra_detail = entry["extra_detail"]
                color = INTERACTION_COLORS.get(interaction_type, "#444444")

                detail_text = f"{interaction_type} {protein_label}"
                if extra_detail:
                    detail_text += f" ({extra_detail})"

                ax.text(
                    0.04,
                    y,
                    "■",
                    color=color,
                    fontsize=11,
                    va="top",
                    transform=ax.transAxes,
                )
                ax.text(
                    0.08,
                    y,
                    detail_text,
                    fontsize=10,
                    va="top",
                    transform=ax.transAxes,
                )
                y -= 0.05

            y -= 0.03

        plt.tight_layout()
        plt.savefig(output_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

    def _merge_main_and_legend(self, main_png, legend_png, output_png):
        main_img = Image.open(main_png).convert("RGB")
        legend_img = Image.open(legend_png).convert("RGB")

        total_width = main_img.width + legend_img.width
        total_height = max(main_img.height, legend_img.height)

        merged = Image.new("RGB", (total_width, total_height), "white")
        merged.paste(main_img, (0, 0))
        merged.paste(legend_img, (main_img.width, 0))
        merged.save(output_png)

        os.remove(main_png)
        os.remove(legend_png)

        return output_png