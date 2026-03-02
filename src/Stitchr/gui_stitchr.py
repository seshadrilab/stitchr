# -*- coding: utf-8 -*-

"""
gui_stitchr.py

A graphical user interface for stitchr, powered by PySide6

"""

import os
import sys
import warnings
import collections as coll

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QTextEdit, QComboBox, QCheckBox, QPushButton,
    QFileDialog, QMessageBox, QScrollArea, QSizePolicy,
)
from PySide6.QtGui import QFont, QShortcut, QKeySequence
from PySide6.QtCore import Qt

from . import stitchrfunctions as fxn
from . import stitchr as st
from . import thimble as th
from . import seqdisplay as sd


__version__ = '1.3.2'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def read_fasta_box(fasta_text):
    """
    :param fasta_text: Contents of text box containing FASTA format text
    """
    header, seq = None, []
    for line in fasta_text:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield header.replace('>', ''), ''.join(seq)
            header, seq = line, []
        else:
            seq.append(line)
    if header:
        yield header, ''.join(seq)


def switch_receptors(text, current_receptor):
    """
    :param text: a str containing a reference to a receptor which needs to be replaced
    :param current_receptor: str of 'TRA/TRB' or 'TRG/TRD', which indicates the required direction of change
    :return: that string with all references to a given locus replaced with the right one
    """
    switched_receptor = text
    changes = [
        ('Alpha', 'Gamma'),
        ('Beta', 'Delta'),
        ('TRA', 'TRG'),
        ('TRB', 'TRD'),
        ('AB', 'GD'),
        ('BA', 'DG')
    ]
    for change in changes:
        if current_receptor == 'TRA/TRB':
            switched_receptor = switched_receptor.replace(change[0], change[1])
        elif current_receptor == 'TRG/TRD':
            switched_receptor = switched_receptor.replace(change[1], change[0])
        else:
            raise ValueError('Unknown receptor configuration detected: ' + current_receptor)
    return switched_receptor


class StitchrWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("stitchr")

        self.extra_gene_text = ">TCRgenename*01\nATG\n"
        self.preferred_button_default = 'Preferred allele file'

        self.linkers = fxn.get_linker_dict()
        self.link_choices = list(self.linkers.keys()) + ['Custom']
        self.species_list = fxn.find_species_covered()
        self.receptor = 'TRA/TRB'
        self.link_orders = {'TRA/TRB': ['AB', 'BA'],
                            'TRG/TRD': ['GD', 'DG']}
        self.examples_path = fxn.gui_examples_dir
        self.preferred = ''
        self.outputs = coll.defaultdict()

        # All named widgets stored here for change_receptors / generic access
        self.widgets = {}

        self._build_ui()
        self._connect_signals()

        QShortcut(QKeySequence(Qt.Key_Escape), self, self.close)

    # ------------------------------------------------------------------ #
    #  UI construction helpers                                             #
    # ------------------------------------------------------------------ #

    def _label(self, text, key=None, font_size=12, font_name='Arial'):
        lbl = QLabel(text)
        lbl.setFont(QFont(font_name, font_size))
        if key:
            self.widgets[key] = lbl
        return lbl

    def _line_edit(self, key, text=''):
        w = QLineEdit(text)
        self.widgets[key] = w
        return w

    def _text_edit(self, key, text='', height_lines=5, mono=True):
        w = QTextEdit()
        if mono:
            w.setFont(QFont('Courier New', 10))
        if text:
            w.setPlainText(text)
        w.setFixedHeight(height_lines * 18)
        self.widgets[key] = w
        return w

    def _combo(self, key, items, default=None):
        w = QComboBox()
        w.addItems(items)
        if default and default in items:
            w.setCurrentText(default)
        self.widgets[key] = w
        return w

    def _checkbox(self, key, text, font_size=12):
        w = QCheckBox(text)
        w.setFont(QFont('Arial', font_size))
        self.widgets[key] = w
        return w

    def _button(self, key, text, font_size=11):
        w = QPushButton(text)
        w.setFont(QFont('Arial', font_size))
        self.widgets[key] = w
        return w

    @staticmethod
    def _row(*widgets):
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        for w in widgets:
            if isinstance(w, QWidget):
                row.addWidget(w)
            else:
                row.addLayout(w)
        return row

    # ------------------------------------------------------------------ #
    #  Full UI build                                                       #
    # ------------------------------------------------------------------ #

    def _build_ui(self):
        # ---- column 1: general controls ----
        col1 = QWidget()
        v1 = QVBoxLayout(col1)
        v1.setAlignment(Qt.AlignTop)

        v1.addLayout(self._row(
            self._button('Example data', 'Example data'),
            self._button('Reset form', 'Reset form'),
        ))

        self._uploaded_tcr_path = ''
        btn_find = self._button('find_tcr_file', 'Find TCR input file')
        btn_upload = self._button('Upload TCR details', 'Upload TCR details')
        v1.addLayout(self._row(btn_find, btn_upload))

        v1.addWidget(self._label('Species', key='species_label'))
        v1.addWidget(self._combo('species_choice', self.species_list, 'HUMAN'))

        v1.addWidget(self._button('change_receptor', 'Change to TRG/TRD'))

        v1.addWidget(self._label('Additional genes', key='additional_genes_label'))
        v1.addWidget(self._text_edit('additional_genes', self.extra_gene_text, height_lines=3))

        btn_preferred = self._button('preferred_allele_button', self.preferred_button_default)
        v1.addWidget(btn_preferred)

        linker_row = QHBoxLayout()
        linker_row.addWidget(self._checkbox('chk_linker', 'Link chains'))
        linker_row.addWidget(self._combo('linker_choice', self.link_choices, self.link_choices[0]))
        v1.addLayout(linker_row)

        v1.addWidget(self._line_edit('custom_linker'))
        self.widgets['custom_linker'].setVisible(False)
        self.widgets['custom_linker'].setPlaceholderText('Custom linker sequence')

        order_row = QHBoxLayout()
        order_row.addWidget(self._label('Link order'))
        order_row.addWidget(self._combo('link_order_choice',
                                        self.link_orders[self.receptor],
                                        self.link_orders[self.receptor][1]))
        v1.addLayout(order_row)

        v1.addWidget(self._checkbox('chk_seamless', 'CDR3 flanking nucleotides (20)'))
        v1.addWidget(self._checkbox('chk_restriction', 'Add Restriction Sites (BamHI, SalI)'))

        run_btn = self._button('Run Stitchr', 'Run Stitchr', font_size=18)
        run_btn.setMinimumHeight(50)
        v1.addWidget(run_btn)

        export_row = QHBoxLayout()
        export_row.addWidget(self._button('Export output', 'Export output'))
        export_row.addWidget(self._button('Exit', 'Exit'))
        v1.addLayout(export_row)

        v1.addWidget(self._label('Linked out', key='linked_out_text'))
        v1.addWidget(self._text_edit('linked_out', height_lines=10))
        v1.addWidget(self._label('Linked log', key='linked_log_text'))
        v1.addWidget(self._text_edit('linked_log', height_lines=5))

        # ---- column 2: Alpha/Gamma chain ----
        col2 = QWidget()
        v2 = QVBoxLayout(col2)
        v2.setAlignment(Qt.AlignTop)

        v2.addWidget(self._label('Alpha chain TCR', key='TR1_title_text', font_size=16))

        v2.addWidget(self._label('TRAV gene name*', key='TR1_V_text'))
        v2.addWidget(self._line_edit('TR1V'))
        v2.addWidget(self._label('TRAJ gene name*', key='TR1_J_text'))
        v2.addWidget(self._line_edit('TR1J'))
        v2.addWidget(self._label('TRA CDR3 junction* (nt/aa)', key='TR1_CDR3_text'))
        v2.addWidget(self._line_edit('TR1_CDR3'))
        v2.addWidget(self._label('TRA arbitrary name', key='TR1_name_text'))
        v2.addWidget(self._line_edit('TR1_name'))

        lc_row = QHBoxLayout()
        lc_row.addWidget(self._label('TRA alternative leader', key='TR1_l_title_text'))
        lc_row.addWidget(self._label('TRAC gene name', key='TR1_c_title_text'))
        v2.addLayout(lc_row)
        lc_inputs = QHBoxLayout()
        lc_inputs.addWidget(self._line_edit('TR1_leader'))
        lc_inputs.addWidget(self._line_edit('TR1C'))
        v2.addLayout(lc_inputs)

        prime_labels = QHBoxLayout()
        prime_labels.addWidget(self._label("5' chain append"))
        prime_labels.addWidget(self._label("3' chain append"))
        v2.addLayout(prime_labels)
        prime_inputs = QHBoxLayout()
        prime_inputs.addWidget(self._line_edit('TR1_5_prime_seq'))
        prime_inputs.addWidget(self._line_edit('TR1_3_prime_seq'))
        v2.addLayout(prime_inputs)

        out1_row = QHBoxLayout()
        out1_row.addWidget(self._label('TRA out', key='TR1_out_title_text'))
        hl1 = self._button('TR1_Highlight', 'Highlight')
        hl1.setEnabled(False)
        out1_row.addWidget(hl1)
        v2.addLayout(out1_row)
        v2.addWidget(self._text_edit('TR1_out', height_lines=20))
        v2.addWidget(self._label('TRA log', key='TR1_log_title_text'))
        v2.addWidget(self._text_edit('TR1_log', height_lines=5))

        # ---- column 3: Beta/Delta chain ----
        col3 = QWidget()
        v3 = QVBoxLayout(col3)
        v3.setAlignment(Qt.AlignTop)

        v3.addWidget(self._label('Beta chain TCR', key='TR2_title_text', font_size=16))

        v3.addWidget(self._label('TRBV gene name*', key='TR2_V_text'))
        v3.addWidget(self._line_edit('TR2V'))
        v3.addWidget(self._label('TRBJ gene name*', key='TR2_J_text'))
        v3.addWidget(self._line_edit('TR2J'))
        v3.addWidget(self._label('TRB CDR3 junction* (nt/aa)', key='TR2_CDR3_text'))
        v3.addWidget(self._line_edit('TR2_CDR3'))
        v3.addWidget(self._label('TRB arbitrary name', key='TR2_name_text'))
        v3.addWidget(self._line_edit('TR2_name'))

        lc2_row = QHBoxLayout()
        lc2_row.addWidget(self._label('TRB alternative leader', key='TR2_l_title_text'))
        lc2_row.addWidget(self._label('TRBC gene name', key='TR2_c_title_text'))
        v3.addLayout(lc2_row)
        lc2_inputs = QHBoxLayout()
        lc2_inputs.addWidget(self._line_edit('TR2_leader'))
        lc2_inputs.addWidget(self._line_edit('TR2C'))
        v3.addLayout(lc2_inputs)

        prime2_labels = QHBoxLayout()
        prime2_labels.addWidget(self._label("5' chain append"))
        prime2_labels.addWidget(self._label("3' chain append"))
        v3.addLayout(prime2_labels)
        prime2_inputs = QHBoxLayout()
        prime2_inputs.addWidget(self._line_edit('TR2_5_prime_seq'))
        prime2_inputs.addWidget(self._line_edit('TR2_3_prime_seq'))
        v3.addLayout(prime2_inputs)

        out2_row = QHBoxLayout()
        out2_row.addWidget(self._label('TRB out', key='TR2_out_title_text'))
        hl2 = self._button('TR2_Highlight', 'Highlight')
        hl2.setEnabled(False)
        out2_row.addWidget(hl2)
        v3.addLayout(out2_row)
        v3.addWidget(self._text_edit('TR2_out', height_lines=20))
        v3.addWidget(self._label('TRB log', key='TR2_log_title_text'))
        v3.addWidget(self._text_edit('TR2_log', height_lines=5))

        # ---- scroll areas for each column ----
        def scrolled(widget):
            sa = QScrollArea()
            sa.setWidgetResizable(True)
            sa.setWidget(widget)
            return sa

        # ---- main layout ----
        central = QWidget()
        main_row = QHBoxLayout(central)
        main_row.addWidget(scrolled(col1))
        main_row.addWidget(scrolled(col2))
        main_row.addWidget(scrolled(col3))
        self.setCentralWidget(central)
        self.resize(1400, 900)

    # ------------------------------------------------------------------ #
    #  Signal connections                                                  #
    # ------------------------------------------------------------------ #

    def _connect_signals(self):
        self.widgets['Example data'].clicked.connect(self._on_example_data)
        self.widgets['Reset form'].clicked.connect(self._on_reset_form)
        self.widgets['find_tcr_file'].clicked.connect(self._on_find_tcr_file)
        self.widgets['Upload TCR details'].clicked.connect(self._on_upload_tcr_details)
        self.widgets['change_receptor'].clicked.connect(self._on_change_receptor)
        self.widgets['preferred_allele_button'].clicked.connect(self._on_preferred_alleles)
        self.widgets['linker_choice'].currentTextChanged.connect(self._on_linker_choice)
        self.widgets['Run Stitchr'].clicked.connect(self._on_run_stitchr)
        self.widgets['Export output'].clicked.connect(self._on_export_output)
        self.widgets['Exit'].clicked.connect(self.close)
        self.widgets['TR1_Highlight'].clicked.connect(self._on_tr1_highlight)
        self.widgets['TR2_Highlight'].clicked.connect(self._on_tr2_highlight)

    # ------------------------------------------------------------------ #
    #  Receptor switching                                                  #
    # ------------------------------------------------------------------ #

    def _change_receptors(self):
        """Updates all labeled fields when swapping between a/b and g/d TCRs."""
        new_receptor = switch_receptors(self.receptor, self.receptor)

        # Update link order combo
        lo = self.widgets['link_order_choice']
        lo.clear()
        lo.addItems(self.link_orders[new_receptor])
        lo.setCurrentText(self.link_orders[new_receptor][1])

        # Update all _text labels and the change_receptor button
        text_keys = [k for k in self.widgets if k.endswith('_text')]
        text_keys.append('change_receptor')

        for key in text_keys:
            w = self.widgets[key]
            if isinstance(w, QLabel):
                current = w.text()
            elif isinstance(w, QPushButton):
                current = w.text()
            else:
                continue

            if key != 'change_receptor':
                new_text = switch_receptors(current, self.receptor)
            else:
                new_text = switch_receptors(current, new_receptor)

            w.setText(new_text)

        self.receptor = new_receptor

    # ------------------------------------------------------------------ #
    #  Value helpers                                                       #
    # ------------------------------------------------------------------ #

    def _val(self, key):
        """Get the current string value of a widget by key."""
        w = self.widgets[key]
        if isinstance(w, QLineEdit):
            return w.text()
        elif isinstance(w, QTextEdit):
            return w.toPlainText()
        elif isinstance(w, QComboBox):
            return w.currentText()
        elif isinstance(w, QCheckBox):
            return w.isChecked()
        return ''

    def _set(self, key, value):
        """Set the value of a widget by key."""
        w = self.widgets[key]
        if isinstance(w, QLineEdit):
            w.setText(str(value) if value is not None else '')
        elif isinstance(w, QTextEdit):
            w.setPlainText(str(value) if value is not None else '')
        elif isinstance(w, QComboBox):
            w.setCurrentText(str(value))
        elif isinstance(w, QCheckBox):
            w.setChecked(bool(value))
        elif isinstance(w, QLabel):
            w.setText(str(value))
        elif isinstance(w, QPushButton):
            w.setText(str(value))

    # ------------------------------------------------------------------ #
    #  Event handlers                                                      #
    # ------------------------------------------------------------------ #

    def _on_example_data(self):
        species = self._val('species_choice')
        example_files = [x for x in os.listdir(self.examples_path)
                         if x.endswith('tsv') and x[0] not in ['.', '~', '_']]
        example_matches = [x for x in example_files
                           if species in x.upper() and self.receptor.replace('/', '-') in x]

        if len(example_matches) == 0:
            QMessageBox.information(self, 'No examples',
                                    'No example ' + self.receptor + ' TCR example files available for species ' + species + '.')
        else:
            if len(example_matches) > 1:
                QMessageBox.information(self, 'Multiple examples',
                                        'More than one ' + self.receptor + 'TCR example files available for species ' + species + ':\n'
                                        'using the first alphabetically.')
                example_matches.sort()

            self._upload_tcr_details(os.path.join(self.examples_path, example_matches[0]), species)

        self.widgets['TR1_Highlight'].setEnabled(False)
        self.widgets['TR2_Highlight'].setEnabled(False)

    def _on_reset_form(self):
        self._set('species_choice', 'HUMAN')

        fields_to_reset = [
            'TR1V', 'TR1J', 'TR1_CDR3', 'TR1_name', 'TR1_leader', 'TR1C',
            'TR1_5_prime_seq', 'TR1_3_prime_seq', 'TR1_out',
            'TR2V', 'TR2J', 'TR2_CDR3', 'TR2_name', 'TR2_leader', 'TR2C',
            'TR2_5_prime_seq', 'TR2_3_prime_seq', 'TR2_out',
        ]
        for field in fields_to_reset:
            self._set(field, '')

        for field in ['linked_out', 'linked_log', 'TR1_log', 'TR2_log']:
            self._set(field, '')

        self.outputs = coll.defaultdict()

        lo = self.widgets['link_order_choice']
        lo.clear()
        lo.addItems(self.link_orders[self.receptor])
        lo.setCurrentText(self.link_orders[self.receptor][1])

        self._set('additional_genes', self.extra_gene_text)
        self._set('preferred_allele_button', self.preferred_button_default)
        self.preferred = ''
        self._uploaded_tcr_path = ''

        self.widgets['TR1_Highlight'].setEnabled(False)
        self.widgets['TR2_Highlight'].setEnabled(False)

    def _on_find_tcr_file(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Find TCR input file')
        if path:
            self._uploaded_tcr_path = path

    def _on_upload_tcr_details(self):
        if not self._uploaded_tcr_path:
            QMessageBox.information(self, 'No file', "Please use 'Find TCR input file' button to try again.")
            return
        self._upload_tcr_details(self._uploaded_tcr_path, self._val('species_choice'))

    def _upload_tcr_details(self, path_to_file, stated_species):
        if not os.path.isfile(path_to_file):
            QMessageBox.information(self, 'TCR file not found',
                                    "Please use 'Find TCR input file' button to try again.")
            return

        species_inference = fxn.infer_species(path_to_file)
        if species_inference:
            inferred_species = species_inference
            self._set('species_choice', inferred_species)
        else:
            inferred_species = stated_species
            QMessageBox.information(self, '', "Cannot infer species name from file name:\nplease set manually.")

        with open(path_to_file, 'r') as in_file:
            line_count = 0
            for line in in_file:
                bits = line.replace('\n', '').replace('\r', '').split('\t')

                if line_count == 0:
                    if bits != th.in_headers[self.receptor]:
                        switched = switch_receptors(self.receptor, self.receptor)
                        if bits == th.in_headers[switched]:
                            self._change_receptors()
                        else:
                            QMessageBox.information(self, 'TCR file error',
                                                    "Input TCR file doesn't have expected columns.\n"
                                                    "Please refer to template and try again.")
                            break

                elif line_count == 1:
                    for x in range(len(th.in_headers[self.receptor])):
                        if x == 0:
                            self._set('TR1_name', bits[x])
                            self._set('TR2_name', bits[x])
                        elif 'Link' not in th.in_headers[self.receptor][x]:
                            self._set(th.locus_to_trx(th.in_headers[self.receptor][x]), bits[x])
                        elif th.in_headers[self.receptor][x] == 'Link_order':
                            self._set('chk_linker', True)
                            if bits[x]:
                                if bits[x] in ('AB', 'GD'):
                                    position = 0
                                elif bits[x] in ('BA', 'DG'):
                                    position = 1
                                else:
                                    warnings.warn("Invalid link order: " + bits[x])
                                    position = 1
                                lo = self.widgets['link_order_choice']
                                lo.clear()
                                lo.addItems(self.link_orders[self.receptor])
                                lo.setCurrentIndex(position)
                        elif th.in_headers[self.receptor][x] == 'Linker':
                            self._set('chk_linker', True)
                            if bits[x]:
                                if bits[x] in self.linkers:
                                    self._set('linker_choice', bits[x])
                                else:
                                    self._set('linker_choice', 'Custom')
                                    self._set('custom_linker', bits[x])
                                    self.widgets['custom_linker'].setVisible(True)

                elif line_count > 1:
                    warnings.warn("More than one data line detected in input TCR file. Ignoring lines after first.")
                    break

                line_count += 1

    def _on_change_receptor(self):
        self._change_receptors()

    def _on_preferred_alleles(self):
        path, _ = QFileDialog.getOpenFileName(self, 'Select preferred allele file')
        if path:
            self.preferred = path
            self._set('preferred_allele_button', os.path.basename(path))

    def _on_linker_choice(self, text):
        self._set('chk_linker', True)
        if text == 'Custom':
            self.widgets['custom_linker'].setVisible(True)
        else:
            self.widgets['custom_linker'].setVisible(False)

    def _tidy_values(self, ref_chain, v):
        for field in ['V', 'J', '_CDR3', '_leader', 'C']:
            k = ref_chain + field
            if v.get(k):
                v[k] = v[k].upper()
        return v

    def _gather_values(self):
        """Collect all current widget values into a dict keyed by widget name."""
        v = {}
        for key, w in self.widgets.items():
            if isinstance(w, QLineEdit):
                v[key] = w.text()
            elif isinstance(w, QTextEdit):
                v[key] = w.toPlainText()
            elif isinstance(w, QComboBox):
                v[key] = w.currentText()
            elif isinstance(w, QCheckBox):
                v[key] = w.isChecked()
        return v

    def _on_run_stitchr(self):
        warning_msgs = coll.defaultdict(str)
        self._set('linked_out', '')
        self._set('linked_log', '')
        self.widgets['Run Stitchr'].setEnabled(False)
        QApplication.processEvents()

        values = self._gather_values()
        values_raw = dict(values)

        codons = fxn.get_optimal_codons('', values['species_choice'])
        self.outputs = coll.defaultdict()

        # Additional genes
        additional_genes_text = values['additional_genes']
        if additional_genes_text != self.extra_gene_text + '\n' and additional_genes_text != self.extra_gene_text:
            self.outputs['additional_fastas_raw'] = [
                x for x in read_fasta_box(additional_genes_text.split('\n') + ['>\n'])
            ][:-1]

            if len(list(set([x[0] for x in self.outputs['additional_fastas_raw']]))) != \
                    len(self.outputs['additional_fastas_raw']):
                self._set('additional_genes',
                          "Multiple FASTAs detected with the same identifier name.\n"
                          "Additional genes ignored; correct and retry")

            self.outputs['additional_fastas'] = []
            for extra_gene in self.outputs['additional_fastas_raw']:
                extra_gene = [x.upper() for x in extra_gene]
                if '*' in extra_gene[0]:
                    self.outputs['additional_fastas'].append(extra_gene)
                else:
                    self.outputs['additional_fastas'].append((extra_gene[0] + '*01', extra_gene[1]))
                if not fxn.dna_check(extra_gene[1]):
                    if fxn.dna_check(extra_gene[1]) != self.extra_gene_text:
                        warnings.warn("Warning: user-provided gene " + extra_gene[0] +
                                      " contains non-DNA sequences.")

        seamless = values['chk_seamless']
        parts = []
        self.outputs['parts'] = parts
        Seq_5 = ""
        Seq_3 = ""
        restriction = False

        if values['chk_restriction']:
            if values['chk_linker']:
                Seq_5 = "GGATCC"
                Seq_3 = "GTCGAC"
            else:
                restriction = True

        convert_chains = {'TRA/TRB': {'TR1': 'TRA', 'TR2': 'TRB'},
                          'TRG/TRD': {'TR1': 'TRG', 'TR2': 'TRD'}}

        for ref_chain in ['TR1', 'TR2']:
            chain = convert_chains[self.receptor][ref_chain]
            self._set(ref_chain + '_out', '')
            self._set(ref_chain + '_log', '')

            with warnings.catch_warnings(record=True) as chain_log:
                warnings.simplefilter("always")

                if values[ref_chain + 'V'] and values[ref_chain + 'J'] and values[ref_chain + '_CDR3']:
                    values = self._tidy_values(ref_chain, values)

                    try:
                        tcr_dat, functionality, partial, frame_dat = fxn.get_imgt_data(
                            chain, st.gene_types, values['species_choice'])

                        if 'additional_fastas' in self.outputs:
                            for extra_gene in self.outputs['additional_fastas']:
                                gene, allele = extra_gene[0].split('*')
                                for gene_type in tcr_dat:
                                    if gene not in tcr_dat[gene_type]:
                                        tcr_dat[gene_type][gene] = coll.defaultdict(list)
                                    if allele in tcr_dat[gene_type][gene]:
                                        raise warnings.warn(
                                            "User provided gene/allele combination " + extra_gene[0] +
                                            " already exists in TCR germline data. Please change and try.")
                                    else:
                                        tcr_dat[gene_type][gene][allele] = extra_gene[1].upper()
                                        functionality[gene][allele] = '?'

                        preferred = fxn.get_preferred_alleles(
                            self.preferred, list(fxn.regions.values()), tcr_dat, partial, chain
                        ) if self.preferred else ''

                        tcr_bits = {
                            'v': values[ref_chain + 'V'],
                            'j': values[ref_chain + 'J'],
                            'cdr3': values[ref_chain + '_CDR3'],
                            'skip_c_checks': False,
                            'species': values['species_choice'],
                            'seamless': seamless,
                            'name': values[ref_chain + '_name'].replace(' ', '_'),
                            'l': values[ref_chain + '_leader'],
                            'c': values[ref_chain + 'C'],
                            '5_prime_seq': values[ref_chain + '_5_prime_seq'],
                            '3_prime_seq': values[ref_chain + '_3_prime_seq'],
                        }

                        if 'additional_fastas' in self.outputs:
                            tcr_bits['skip_c_checks'] = True

                        tcr_bits = fxn.autofill_input(tcr_bits, chain)

                        mouse_c = ''
                        if values['species_choice'] == 'HUMAN' and not values_raw[ref_chain + 'C']:
                            if chain == 'TRA':
                                mouse_c = ('TRAC*00', 'gacattcagaacccggaaccggctgtataccagctgaaggacccccgatctcaggatagtactctgtgcctgttcaccgactttgatagtcagatcaatgtgcctaaaaccatggaatccggaacttttattaccgacaagtgcgtgctggatatgaaagccatggacagtaagtcaaacggcgccatcgcttggagcaatcagacatccttcacttgccaggatatcttcaaggagaccaacgcaacatacccatcctctgacgtgccctgtgatgccaccctgacagagaagtctttcgaaacagacatgaacctgaattttcagaatctgagcgtgatgggcctgagaatcctgctgctgaaggtcgctgggtttaatctgctgatgacactgcggctgtggtcctca'.upper())
                            elif chain == 'TRB':
                                mouse_c = ('TRBC*00', 'gaagatctacgtaacgtgacaccacccaaagtctcactgtttgagcctagcaaggcagaaattgccaacaagcagaaggccaccctggtgtgcctggcaagagggttctttccagatcacgtggagctgtcctggtgggtcaacggcaaagaagtgcattctggggtctgcaccgacccccaggcttacaaggagagtaattactcatattgtctgtcaagccggctgagagtgtccgccacattctggcacaaccctaggaatcatttccgctgccaggtccagtttcacggcctgagtgaggaagataaatggccagaggggtcacctaagccagtgacacagaacatcagcgcagaagcctggggacgagcagactgtggcattactagcgcctcctatcatcagggcgtgctgagcgccactatcctgtacgagattctgctgggaaaggccaccctgtatgctgtgctggtctccggcctggtgctgatggccatggtcaagaaaaagaactct'.upper())

                        self.outputs[ref_chain + '_out_list'], \
                        self.outputs[ref_chain + '_stitched'], \
                        self.outputs[ref_chain + '_offset'], region, check = st.stitch(
                            tcr_bits, tcr_dat, functionality, partial, codons, 3, preferred, mouse_c, frame_dat, restriction)

                        self.outputs[ref_chain + '_out_str'] = '|'.join(self.outputs[ref_chain + '_out_list'])
                        self.outputs[ref_chain + '_fasta'] = fxn.fastafy(
                            'nt|' + self.outputs[ref_chain + '_out_str'],
                            self.outputs[ref_chain + '_stitched'])
                        self._set(ref_chain + '_out', self.outputs[ref_chain + '_fasta'])

                        if not values['chk_linker']:
                            seq = fxn.translate_nt(self.outputs[ref_chain + '_stitched'])
                            check_seq = seq[len(seq) - 3] if restriction else seq[len(seq) - 1]
                            if check_seq != '*':
                                check.append("Warning: Stop codon expected, but not found at end of sequence.")
                            else:
                                check.append("Check: Stop codon successfully located.")

                        parts.append(region)
                        warning_msgs[ref_chain + '_out'] = '\n'.join([str(check[x]) for x in range(len(check))])

                    except Exception as message:
                        warning_msgs[ref_chain + '_out'] = str(message)

                elif values[ref_chain + 'V'] or values[ref_chain + 'J'] or values[ref_chain + '_CDR3']:
                    warnings.warn('V gene, J gene, and CDR3 sequence are all required to stitch a TCR chain.')

            warning_msgs[ref_chain + '_out'] += '\n'.join([
                str(chain_log[x].message) for x in range(len(chain_log))
                if 'DeprecationWarning' not in str(chain_log[x].category)
            ])

            if not values['chk_linker']:
                self.widgets[ref_chain + '_Highlight'].setEnabled(True)
            self._set(ref_chain + '_log', warning_msgs[ref_chain + '_out'])

        # Link chains if requested
        if values['chk_linker']:
            with warnings.catch_warnings(record=True) as link_log:
                warnings.simplefilter("always")
                try:
                    if 'TR1_out_str' in self.outputs and 'TR2_out_str' in self.outputs:
                        link_order = values['link_order_choice']
                        if link_order in ('BA', 'DG'):
                            tr1, tr2 = '2', '1'
                        elif link_order in ('AB', 'GD'):
                            tr1, tr2 = '1', '2'
                        else:
                            raise warnings.warn("Undetermined link order.")

                        self.outputs['linker'] = values['linker_choice']
                        if self.outputs['linker'] == 'Custom':
                            custom = values['custom_linker']
                            if custom:
                                self.linkers['Custom'] = custom
                            else:
                                self._set('linked_log',
                                          "Cannot output linked sequence: custom linker chosen, but not provided")

                        self.outputs['linker_seq'] = fxn.get_linker_seq(self.outputs['linker'], self.linkers)
                        self.outputs['linked'] = (Seq_5 +
                                                  self.outputs['TR' + tr1 + '_stitched'] +
                                                  self.outputs['linker_seq'] +
                                                  self.outputs['TR' + tr2 + '_stitched'] +
                                                  Seq_3)
                        self.outputs['linked_header'] = '_'.join([
                            self.outputs['TR' + tr1 + '_out_str'],
                            self.outputs['linker'],
                            self.outputs['TR' + tr2 + '_out_str'],
                        ])
                        self.outputs['linked_fasta'] = fxn.fastafy(
                            self.outputs['linked_header'], self.outputs['linked'])
                        self._set('linked_out', self.outputs['linked_fasta'])

                        seq = fxn.translate_nt(self.outputs['linked'])
                        check_seq = seq[len(seq) - 3] if values['chk_restriction'] else seq[len(seq) - 1]
                        if check_seq != '*':
                            warning_msgs['linked_out'] += "Warning: Stop codon expected, but not found at end of sequence."
                        else:
                            warning_msgs['linked_out'] += "Check: Stop codon successfully located."

                    else:
                        raise warnings.warn("Two valid chains required for linking.")

                except Exception as message:
                    warning_msgs['linked_out'] += str(message)

            warning_msgs['linked_out'] += ''.join([
                str(link_log[x].message) for x in range(len(link_log))
                if 'DeprecationWarning' not in str(link_log[x].category)
            ])

            if warning_msgs['linked_out']:
                self._set('linked_log', warning_msgs['linked_out'])

            if not seamless and 'linked' in self.outputs:
                sd.display(self.outputs['linked'], parts, fxn.translate_nt(self.outputs['linker_seq']), True)

        self.widgets['Run Stitchr'].setEnabled(True)

    def _on_export_output(self):
        if not self.outputs:
            return
        out_str = ''
        for chain in ['TR1', 'TR2']:
            if chain + '_fasta' in self.outputs:
                out_str += self.outputs[chain + '_fasta']

        values = self._gather_values()
        if values['chk_linker'] and ('TR1_out_str' in self.outputs and 'TR2_out_str' in self.outputs):
            out_str += self.outputs['linked_fasta']

        if not out_str:
            return

        out_file, _ = QFileDialog.getSaveFileName(
            self, 'Export output', '', 'FASTA files (*.fasta)')
        if out_file:
            with open(out_file, 'w') as f:
                f.write(out_str)

    def _on_tr1_highlight(self):
        if 'TR1_stitched' in self.outputs:
            sd.display(self.outputs['TR1_stitched'], self.outputs.get('parts', []), "")
        self.widgets['TR1_Highlight'].setEnabled(False)

    def _on_tr2_highlight(self):
        if 'TR2_stitched' in self.outputs:
            sd.display(self.outputs['TR2_stitched'], self.outputs.get('parts', []), "")
        self.widgets['TR2_Highlight'].setEnabled(False)


def main():
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    window = StitchrWindow()
    window.show()
    sys.exit(app.exec())
