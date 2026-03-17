import re
from . import stitchrfunctions as fxn
from PySide6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QTextEdit, QPushButton, QApplication
from PySide6.QtGui import QTextCharFormat, QColor, QTextCursor, QFont
from PySide6.QtCore import Qt


def get_partlist(chains, linker):
    """
    param chains: A list of dictionary of amino acids and their associated lables (v, j, c, CDR3)
    param linker: A string of amino acids that connects the two TCRs in sequence
    return: names, a list of different parts to highlight.  parts, an ordered list of different amino acid chains associated with names
    """
    TR = 0
    parts = []
    names = []
    for chain in chains:
        TR+=1
        for query in chain:
            names.append(str(TR) + "_" + query)
            parts.append(chain[query])
    names.append("Linker")
    parts.append(linker)
    names.append("Start")
    parts.append('M')
    names.append("End")
    parts.append('*')
    return names, parts


def get_indexes(seq, name, part, a=1):
    """
    param seq: A string amino acid sequence representing our stitched TCR
    param name: A list of named parts to highlight
    param part: A list of different amino acid sequences that correspond to names in name list
    param a: defaults to 1.  Is used as a modifier depending on whether in NT mode (a=3) or AA mode (a=1), changes positions to highlight
    return: zipped list of named sequences that show up in the stitched TCR and the positions to highight (start and end indexes)
    """
    parts = list(zip(name, part))
    names = []
    index1s = []
    index2s = []
    for name, part in parts:
        if part in seq:
            if name == 'Start' or name == 'End':
                start = seq.index(part)
                index1s.append(a * start)
                index2s.append(a * (start + 1))
                names.append(name)
            else:
                for m in re.finditer(part, seq):
                    index1s.append(a * m.start())
                    index2s.append(a * m.end())
                    names.append(name)
    indexes = list(zip(names, index1s, index2s))
    return indexes


# Map region name patterns to (foreground, background) colors
_COLOR_MAP = [
    ("_cdr3",   ("white",  "black")),
    ("_l",      ("white",  "purple")),
    ("Linker",  ("white",  "teal")),
    ("Start",   ("white",  "green")),
    ("End",     ("white",  "red")),
    ("_v",      ("black",  "orange")),
    ("_c",      ("black",  "pink")),
    ("_j",      ("white",  "brown")),
]


def _make_fmt(fg, bg, bold=False):
    fmt = QTextCharFormat()
    fmt.setForeground(QColor(fg))
    fmt.setBackground(QColor(bg))
    if bold:
        fmt.setFontWeight(QFont.Bold)
    return fmt


def _get_color(name):
    for pattern, colors in _COLOR_MAP:
        if pattern in name:
            return colors
    return None


def get_highlights(widget, indexes, bold=False):
    """
    param widget: a QTextEdit widget to apply highlights to
    param indexes: a zipped list of (name, start_pos, end_pos) tuples
    param bold: whether to use bold font for highlights
    return: list of QTextEdit.ExtraSelection objects to apply
    """
    selections = []
    for name, start, end in indexes:
        colors = _get_color(name)
        if colors is None:
            continue
        fmt = _make_fmt(colors[0], colors[1], bold)
        sel = QTextEdit.ExtraSelection()
        sel.format = fmt
        cursor = QTextCursor(widget.document())
        cursor.setPosition(start)
        cursor.setPosition(end, QTextCursor.KeepAnchor)
        sel.cursor = cursor
        selections.append(sel)
    return selections


class SeqDisplayDialog(QDialog):

    def __init__(self, nt, parts="", linker="", linked=False, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Sequence Display")

        aa = fxn.translate_nt(nt)
        self._nt = nt
        self._aa = aa

        m_name, m_part = get_partlist(parts, linker)
        self._nt_indexes = get_indexes(aa, m_name, m_part, 3)
        self._aa_indexes = get_indexes(aa, m_name, m_part, 1)
        self._m_indexes = self._aa_indexes

        legend = "leader sequence | Linker sequence | cdr3 sequence | v region | j region | c region | Start | End"
        l_name = ['_' + i for i in legend.split(' | ')]
        l_part = legend.split(' | ')
        self._l_indexes = get_indexes(legend, l_name, l_part)

        font = QFont("Courier New", 10)

        height = 20 if linked else 10

        self._seq_box = QTextEdit()
        self._seq_box.setFont(font)
        self._seq_box.setPlainText(aa)
        self._seq_box.setReadOnly(True)
        self._seq_box.setMinimumWidth(700)
        self._seq_box.setFixedHeight(height * 18)

        self._legend_box = QTextEdit()
        self._legend_box.setFont(font)
        self._legend_box.setPlainText(legend)
        self._legend_box.setReadOnly(True)
        self._legend_box.setFixedHeight(30)

        self._btn_highlight = QPushButton("Highlight")
        self._btn_exit = QPushButton("Exit")
        self._btn_nt = QPushButton("NT")
        self._btn_aa = QPushButton("AA")
        self._btn_aa.setEnabled(False)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_row.addWidget(self._btn_highlight)
        btn_row.addWidget(self._btn_exit)
        btn_row.addWidget(self._btn_nt)
        btn_row.addWidget(self._btn_aa)

        layout = QVBoxLayout()
        layout.addWidget(self._seq_box)
        layout.addWidget(self._legend_box)
        layout.addLayout(btn_row)
        self.setLayout(layout)

        self._btn_highlight.clicked.connect(self._on_highlight)
        self._btn_exit.clicked.connect(self.accept)
        self._btn_nt.clicked.connect(self._on_nt)
        self._btn_aa.clicked.connect(self._on_aa)

    def _on_highlight(self):
        bold_font = QFont("Courier New", 10)
        bold_font.setBold(True)
        seq_sels = get_highlights(self._seq_box, self._m_indexes, bold=True)
        leg_sels = get_highlights(self._legend_box, self._l_indexes, bold=True)
        self._seq_box.setExtraSelections(seq_sels)
        self._legend_box.setExtraSelections(leg_sels)
        self._btn_highlight.setEnabled(False)

    def _on_nt(self):
        self._seq_box.setPlainText(self._nt)
        self._seq_box.setExtraSelections([])
        self._legend_box.setExtraSelections([])
        self._m_indexes = self._nt_indexes
        self._btn_nt.setEnabled(False)
        self._btn_aa.setEnabled(True)
        self._btn_highlight.setEnabled(True)

    def _on_aa(self):
        self._seq_box.setPlainText(self._aa)
        self._seq_box.setExtraSelections([])
        self._legend_box.setExtraSelections([])
        self._m_indexes = self._aa_indexes
        self._btn_nt.setEnabled(True)
        self._btn_aa.setEnabled(False)
        self._btn_highlight.setEnabled(True)


def display(nt, parts="", linker="", linked=False):
    """
    Param nt: A string nt DNA sequence
    Param parts: a dictionary of gene regions and their DNA Amino Acid sequence
    Display: A GUI display of the DNA sequence that highlights different regions
    """
    app = QApplication.instance()
    created_app = app is None
    if created_app:
        import sys
        app = QApplication(sys.argv)

    dialog = SeqDisplayDialog(nt, parts, linker, linked)
    dialog.exec()


def main():
    print("Please use the appropriate 'stitchr', 'thimble', 'gui_stitchr' or 'stitchrdl' command.")
