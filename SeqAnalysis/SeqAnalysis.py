#!/usr/bin/env python3

"""
SeqAnalysis.py

Suite of functions for analyzing phage display data (*.fasta or *.txt/seq) by sorting,
converting to *.fasta (if input file is *.txt/seq), trimming, and translating input sequences.
"""

from tools.path_tools import format_path, strip_filetype, fetch_all_files, create_dir_move_files
from tools.seq_manip import export_batch_fasta
from tools.misc import add_leading_zeros
import re
from pathlib import Path
import os
import xlsxwriter
import json


class PdasSeqAnalysis:
    def __init__(self, input_path, prot_ref, dna_ref, lib_file):
        self.master_path = format_path('', input_path)
        self.prot_ref = prot_ref
        self.dna_ref = dna_ref
        self.lib_file = lib_file
        self.nt_regex = r'[ATCGatcg]+'
        self.trim_motif = r'AAAATG'
        self.trim_length = 231

    def read_seq_folder(self):
        """
        Read in sequencing data (*.fasta, *.seq, *.txt) from a folder and create a batch
        fasta file containing files with sequences.
        """
        formatted_path = format_path('', self.master_path)
        file_types = ['*.seq', '*.txt', '*.fasta']
        files = fetch_all_files(formatted_path, file_types)

        fasta_dict = {}
        for file in files:
            if Path(file).name.endswith('_seq.fasta'):
                continue
            try:
                file_name = strip_filetype(file)
                fasta_name = '>' + add_leading_zeros(file_name, 2)
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"The file at '{formatted_path + file}' was not found.")
            try:
                with open(formatted_path + file, 'r', encoding='utf-8') as f:
                    seq = f.read()
                    seq_clean = seq.replace('\n', '')
                    if not seq_clean:
                        raise ValueError(f"The file '{file}' is empty.")
                    elif re.search(self.nt_regex, seq_clean):
                        fasta_dict[fasta_name] = seq_clean
            except UnicodeDecodeError:
                raise Exception(
                    "There was an error decoding the file. Please check the file encoding.")
            except Exception as e:
                raise Exception(f"An unexpected error occurred: '{e}'")

        return fasta_dict

    def trim_seq(self, input_dict):
        """
        Trim seqs beginning with a desired motif and ending at a desired length.

        The last three characters in the input motif (typically the start codon, ATG) are
        included in output seqs. seqs that don't meet trim requirements are excluded
        from the output.
        """
        trimmed_seqs = {}
        for id, seq in input_dict.items():
            match = re.search(self.trim_motif, seq)
            if match:
                start = match.end() - 3
                trimmed_seq = seq[start:start + self.trim_length + 3]
                if len(trimmed_seq) == self.trim_length + 3:
                    trimmed_seqs[id] = trimmed_seq
        return trimmed_seqs

    def translate_dna_dict_to_protein_dict(self, dna_dict):
        codon_table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'
        }

        protein_dict = {}
        for seq_id, dna_sequence in dna_dict.items():
            protein_sequence = ""
            for i in range(0, len(dna_sequence), 3):
                codon = dna_sequence[i:i+3]
                if codon in codon_table:
                    protein_sequence += codon_table[codon]
                else:
                    protein_sequence += 'X'
            protein_dict[seq_id] = protein_sequence
        return protein_dict

    def count_unique_seqs_with_ids(self, input_dict):
        """
        Find unique sequences, count their occurrences, and relate them to their IDs.
        """
        unique_seq = {}
        for seq_id, sequence in input_dict.items():
            if sequence in unique_seq:
                unique_seq[sequence]['ids'].append(seq_id)
                unique_seq[sequence]['freq'] += 1
            else:
                unique_seq[sequence] = {'ids': [seq_id], 'freq': 1}
        unique_seq = dict(
            sorted(unique_seq.items(), key=lambda item: item[1]['freq'], reverse=True))
        return unique_seq

    def conserved_dict(self, reference_seq, query_dict):
        """
        Replace characters in a dictionary of sequences (DNA or protein) that are conserved in
        the reference sequence with dashes.
        """
        def compare_strings(ref, query):
            diff = ''
            for x, y in zip(ref, query):
                if x == y:
                    diff += '-'
                else:
                    diff += y
            return diff

        new_conserved_dict = {}
        for query_seq, query_info in query_dict.items():
            conserved_seq = compare_strings(reference_seq, query_seq)
            new_conserved_dict[conserved_seq] = query_info
        return new_conserved_dict

    def read_lib(self):
        """
        Read phage display library data (*.json).
        """
        script_path = os.path.abspath(__file__)
        script_dir = os.path.dirname(script_path) + r'\libraries'
        formatted_path = format_path('', script_dir)
        with open(formatted_path + self.lib_file, 'r') as file:
            data = json.load(file)
        return data

    def pd_excel_export(self, nt_dict, aa_dict, nt_dict_conserved, aa_dict_conserved, lib_data):
        """
        Export phage display data in protein and DNA formats as Excel worksheets.
        """
        workbook = xlsxwriter.Workbook(self.master_path + 'results.xlsx')

        # Setup Excel formats
        title_fmt = workbook.add_format({
            'bold': True,
            'font_size': 11,
            'align': 'center',
            'valign': 'vcenter',
            'font_name': 'Segoe UI',
            'font_color': 'white',
            'bg_color': '#5b5b5b',
            'border': 1,
            'border_color': '#999999'
        })
        num_fmt = workbook.add_format({
            'font_size': 7,
            'align': 'center',
            'valign': 'vcenter',
            'font_name': 'Lucida Console',
            'bg_color': '#eeeeee'
        })
        res_fmt = workbook.add_format({
            'font_size': 7,
            'align': 'center',
            'valign': 'vcenter',
            'font_name': 'Lucida Console'
        })
        seq_fmt = workbook.add_format({
            'font_size': 8,
            'align': 'center',
            'valign': 'vcenter',
            'font_name': 'Lucida Console'
        })
        freq_fmt = workbook.add_format({
            'font_size': 8,
            'align': 'center',
            'valign': 'vcenter',
            'font_name': 'Segoe UI'
        })
        id_fmt = workbook.add_format({
            'font_size': 8,
            'align': 'left',
            'valign': 'vcenter',
            'font_name': 'Segoe UI'
        })
        r1_fmt = workbook.add_format({'bg_color': '#BD7191'})
        r2_fmt = workbook.add_format({'bg_color': '#4e78a0'})
        r3_fmt = workbook.add_format({'bg_color': '#f7c654'})

            # Conserved protein sequences
        worksheet1 = workbook.add_worksheet('Protein Seq (Conserved)')
        worksheet1.hide_gridlines(2)
        worksheet1.freeze_panes(4, 0)

        res_len = len(next(iter(aa_dict_conserved))) + 1
        res_col = 0
        res_row = 3
        for i in range(1, res_len):
            worksheet1.write(res_row, res_col, i, num_fmt)
            res_col += 1

        seq_len = len(next(iter(aa_dict_conserved))) - 1
        worksheet1.merge_range(0, 0, 2, seq_len, 'Sequence', title_fmt)
        worksheet1.set_column(0, seq_len, 2.5)
        seq_col = 0
        seq_row = res_row + 1
        for seq in aa_dict_conserved.keys():
            seq_list = list(seq)
            for res in seq_list:
                worksheet1.write(seq_row, seq_col, res, seq_fmt)
                seq_col += 1
            seq_row += 1
            seq_col = 0

        # for lib in lib_data['libraries']:
        #     for region, residues in lib['diversified_residues'].items():
        #             # worksheet1.conditional_format(seq_row, , len(aa_dict_conserved), {'type': 'no_blanks', 'format': r1_fmt})

        freq_col = seq_len + 1
        worksheet1.merge_range(0, freq_col, 2, freq_col, 'Frequency', title_fmt)
        worksheet1.set_column(freq_col, freq_col, 12)
        worksheet1.write(0, freq_col, 'Frequency', title_fmt)
        freq_row = 4
        for seq, data in aa_dict_conserved.items():
            worksheet1.write(freq_row, freq_col, data['freq'], freq_fmt)
            freq_row += 1

        id_col = freq_col + 1
        worksheet1.merge_range(0, id_col, 2, id_col, 'ID', title_fmt)
        worksheet1.set_column(id_col, id_col, 50)
        id_row = 4
        for seq, data in aa_dict_conserved.items():
            ids_str = ', '.join(map(str, data['ids']))
            worksheet1.write(id_row, id_col, ids_str, id_fmt)
            id_row += 1

        # Full protein sequences
        worksheet2 = workbook.add_worksheet('Protein Seq (Full)')
        worksheet2.hide_gridlines(2)
        worksheet2.freeze_panes(4, 0)

        res_len = len(next(iter(aa_dict))) + 1
        res_col = 0
        res_row = 3
        for i in range(1, res_len):
            worksheet2.write(res_row, res_col, i, num_fmt)
            res_col += 1

        seq_len = len(next(iter(aa_dict))) - 1
        worksheet2.merge_range(0, 0, 2, seq_len, 'Sequence', title_fmt)
        worksheet2.set_column(0, seq_len, 2.5)
        seq_col = 0
        seq_row = res_row + 1
        for seq in aa_dict.keys():
            seq_list = list(seq)
            for res in seq_list:
                worksheet2.write(seq_row, seq_col, res, seq_fmt)
                seq_col += 1
            seq_row += 1
            seq_col = 0

        freq_col = seq_len + 1
        worksheet2.merge_range(0, freq_col, 2, freq_col, 'Frequency', title_fmt)
        worksheet2.set_column(freq_col, freq_col, 12)
        worksheet2.write(0, freq_col, 'Frequency', title_fmt)
        freq_row = 4
        for seq, data in aa_dict.items():
            worksheet2.write(freq_row, freq_col, data['freq'], freq_fmt)
            freq_row += 1

        id_col = freq_col + 1
        worksheet2.merge_range(0, id_col, 2, id_col, 'ID', title_fmt)
        worksheet2.set_column(id_col, id_col, 50)
        id_row = 4
        for seq, data in aa_dict.items():
            ids_str = ', '.join(map(str, data['ids']))
            worksheet2.write(id_row, id_col, ids_str, id_fmt)
            id_row += 1

        # Conserved DNA sequences
        worksheet3 = workbook.add_worksheet('DNA Seq (Conserved)')
        worksheet3.hide_gridlines(2)
        worksheet3.freeze_panes(4, 0)

        res_len = len(next(iter(nt_dict_conserved))) + 1
        res_col = 0
        res_row = 3
        for i in range(1, res_len):
            worksheet3.write(res_row, res_col, i, num_fmt)
            res_col += 1

        seq_len = len(next(iter(nt_dict_conserved))) - 1
        worksheet3.merge_range(0, 0, 2, seq_len, 'Sequence', title_fmt)
        worksheet3.set_column(0, seq_len, 2.5)
        seq_col = 0
        seq_row = res_row + 1
        for seq in nt_dict_conserved.keys():
            seq_list = list(seq)
            for res in seq_list:
                worksheet3.write(seq_row, seq_col, res, seq_fmt)
                seq_col += 1
            seq_row += 1
            seq_col = 0

        freq_col = seq_len + 1
        worksheet3.merge_range(0, freq_col, 2, freq_col, 'Frequency', title_fmt)
        worksheet3.set_column(freq_col, freq_col, 12)
        worksheet3.write(0, freq_col, 'Frequency', title_fmt)
        freq_row = 4
        for seq, data in nt_dict_conserved.items():
            worksheet3.write(freq_row, freq_col, data['freq'], freq_fmt)
            freq_row += 1

        id_col = freq_col + 1
        worksheet3.merge_range(0, id_col, 2, id_col, 'ID', title_fmt)
        worksheet3.set_column(id_col, id_col, 50)
        id_row = 4
        for seq, data in nt_dict_conserved.items():
            ids_str = ', '.join(map(str, data['ids']))
            worksheet3.write(id_row, id_col, ids_str, id_fmt)
            id_row += 1

        # Full DNA sequences
        worksheet4 = workbook.add_worksheet('DNA Seq (Full)')
        worksheet4.hide_gridlines(2)
        worksheet4.freeze_panes(4, 0)

        res_len = len(next(iter(nt_dict))) + 1
        res_col = 0
        res_row = 3
        for i in range(1, res_len):
            worksheet4.write(res_row, res_col, i, num_fmt)
            res_col += 1

        seq_len = len(next(iter(nt_dict))) - 1
        worksheet4.merge_range(0, 0, 2, seq_len, 'Sequence', title_fmt)
        worksheet4.set_column(0, seq_len, 2.5)
        seq_col = 0
        seq_row = res_row + 1
        for seq in nt_dict.keys():
            seq_list = list(seq)
            for res in seq_list:
                worksheet4.write(seq_row, seq_col, res, seq_fmt)
                seq_col += 1
            seq_row += 1
            seq_col = 0

        freq_col = seq_len + 1
        worksheet4.merge_range(0, freq_col, 2, freq_col, 'Frequency', title_fmt)
        worksheet4.set_column(freq_col, freq_col, 12)
        worksheet4.write(0, freq_col, 'Frequency', title_fmt)
        freq_row = 4
        for seq, data in nt_dict.items():
            worksheet4.write(freq_row, freq_col, data['freq'], freq_fmt)
            freq_row += 1

        id_col = freq_col + 1
        worksheet4.merge_range(0, id_col, 2, id_col, 'ID', title_fmt)
        worksheet4.set_column(id_col, id_col, 50)
        id_row = 4
        for seq, data in nt_dict.items():
            ids_str = ', '.join(map(str, data['ids']))
            worksheet4.write(id_row, id_col, ids_str, id_fmt)
            id_row += 1

        workbook.close()

    def run_analysis(self, output_path):
        """
        Run the sequence analysis process.
        """
        dna_sequences = self.read_seq_folder()
        raw_ext = ['.seq', '.txt', '.ab1']
        # create_dir_move_files(os.path.join(output_path, 'raw_data'),
        #                       source_dir=output_path,
        #                       filetype_list=raw_ext
        #                       )
        
        trimmed_sequences = self.trim_seq(dna_sequences)
        # for id, seq in trimmed_sequences.items():
        #     export_batch_fasta(id, seq, 'dna_seq', output_path)
            
        protein_sequences = self.translate_dna_dict_to_protein_dict(trimmed_sequences)
        # for id, seq in protein_sequences.items():
        #     export_batch_fasta(id, seq, 'prot_seq', output_path)
            
        unique_dna_sequences = self.count_unique_seqs_with_ids(trimmed_sequences)
        unique_protein_sequences = self.count_unique_seqs_with_ids(protein_sequences)
        
        dna_conserved = self.conserved_dict(self.dna_ref, unique_dna_sequences)
        protein_conserved = self.conserved_dict(self.prot_ref, unique_protein_sequences)
        
        library_data = self.read_lib()
        self.pd_excel_export(unique_dna_sequences,
                             unique_protein_sequences,
                             dna_conserved,
                             protein_conserved,
                             library_data
                             )


if __name__ == "__main__":
    
    analysis = PdasSeqAnalysis(
        input_path=r'C:\Users\spyro\Documents\Work\programming\phage_display_analysis_suite_pdas\SeqAnalysis\files\input',
        prot_ref='MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGGGG',
        dna_ref='ATGCAGATTTTCGTGAAAACCCTTACGGGGAAGACCATCACCCTCGAGGTTGAACCCTCGGATACGATAGAAAATGTAAAGGCCAAGATCCAGGATAAGGAAGGAATTCCTCCTGATCAGCAGAGACTGATCTTTGCTGGCAAGCAGCTGGAAGATGGACGTACTTTGTCTGACTACAATATTCAAAAGGAGTCTACTCTTCATCTTGTGTTGAGACTTCGTGGTGGTGCTAAGAAA',
        lib_file='ub_lib.json'
    )
    analysis.run_analysis(analysis.master_path)
