<h1 align="center">Phage Display Analysis Suite (PDAS)</h1>

<p align="left">Suite of codes for analyzing bulk phage display sequencing data (*.fasta or *.txt/seq) by reading input DNA sequence files, converting to *.fasta (if input file is *.txt/seq), trimming at a specified motif with a desired length, translating DNA sequences to protein sequences, and colour codes diversified regions using built-in or custom formatting.</p>

<p>This program works very simply; locate the folder with your input files and the library file with the desired formatting, The file will be previewed below once selected and the combo box will be populated with the names of all libraries contained within the library file. Select the desired library name, then run.</p>

<p>Please note the text formatting of the library file (*.json), this formatting must be adhered to for custom libraries to work.</p>

## Links

- [Repo](https://github.com/oroscow/phage_display_analysis_suite_pdas "PDAS Repo")

- [Bugs](https://github.com/oroscow/phage_display_analysis_suite_pdas/issues "Issues Page")

## Screenshots

![PDAS GUI](/files/images/screenshot1.png "PDAS GUI")

*Main interface.*

---

![PDAS input files](/files/images/screenshot2.png "PDAS input files")

*Example input files.*

---

![PDAS output files](/files/images/screenshot3.png "PDAS output files")

*Example output file.*

---

![PDAS output folder](/files/images/screenshot4.png "PDAS output folder")

*File organization post-analysis, where original input files remained untouched in the `raw_data` folder. `dna_seq.fasta` and `prot_seq.fasta` are batch fasta files containing the trimmed and translated version of all sequences included in analysis.*

## Built With

- Python
- QML (PyQt6)
- HTML
- Bash

## Future Updates

- [ ] Set default paths for folder/file browsers (first try didn't work).
- [ ] Use a more descriptive name for final output.
- [ ] Incorporate PyRosetta/AlphaFold/etc.

## Author

**Olivia Roscow**

- [Profile](https://github.com/oroscow "Olivia Roscow")
- [Email](mailto:spyrolivia@gmail.com "spyrolivia@gmail.com")
- [Portfolio](https://oroscow.github.io/ "Portfolio")

## Support

Contributions, issues, and feature requests are welcome!
