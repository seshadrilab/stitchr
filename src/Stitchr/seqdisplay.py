import PySimpleGUI as sg
import re
from . import stitchrfunctions as fxn

def get_partlist(chains, linker):
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
    names.append('End')
    parts.append('*')
    return names, parts


def get_indexes(seq, name, part, a=1):
    """
    Input: A string DNA sequence and a dictionary of gene regions and their DNA sequence
    Output: A list of tuples where each gene region occurs in DNA sequence
    """
    parts = list(zip(name, part))
    names = []
    index1s = []
    index2s = []
    for name, part in parts:
        if part in seq:
            if name == 'Start':
                start = seq.index(part)
                index1s.append(f'1.{a*start}')
                index2s.append(f'1.{a*(start + 1)}')
                names.append(name)
            else:
                for m in re.finditer(part, seq):
                    index1s.append(f'1.{a*m.start()}')
                    index2s.append(f'1.{a*m.end()}')
                    names.append(name)
    indexes = list(zip(names, index1s, index2s))
    return indexes


def get_highlights(widget, indexes, fonts):
    for name, index1, index2 in indexes:
                    # Make it so that a value is stored in indexes that will check against defined (v, j, cdr3, c, start, stop)
                    if "_cdr3" in name:
                        widget.tag_config('BLACK', foreground='white', background='black', font=fonts)
                        widget.tag_add('BLACK', index1, index2)
                    elif "_l" in name:
                        widget.tag_config('PURPLE', foreground='white', background='purple', font=fonts)
                        widget.tag_add('PURPLE', index1, index2)
                    elif "Linker" in name:
                        widget.tag_config('BLUE', foreground='white', background='teal', font=fonts)
                        widget.tag_add('BLUE', index1, index2)
                    elif "Start" in name:
                        widget.tag_config('GREEN', foreground='white', background='green', font=fonts)
                        widget.tag_add('GREEN', index1, index2)
                    elif "End" in name:
                        widget.tag_config('RED', foreground='white', background='red', font=fonts)
                        widget.tag_add('RED', index1, index2)
                    elif "_v" in name:
                        widget.tag_config('ORANGE', foreground='black', background='orange', font=fonts)
                        widget.tag_add('ORANGE', index1, index2)
                    elif "_c" in name:
                        widget.tag_config('PINK', foreground='black', background='pink', font=fonts)
                        widget.tag_add('PINK', index1, index2)
                    elif "_j" in name:
                        widget.tag_config('BROWN', foreground='white', background='brown', font=fonts)
                        widget.tag_add('BROWN', index1, index2)
    return widget

def display(nt, parts, linker=[]):
    """
    Input: A string DNA sequence and a dictionary of gene regions and their DNA sequence
    Output: A GUI display of the DNA sequence that highlights different regions
    """
    aa = fxn.translate_nt(nt)
    m_name, m_part = get_partlist(parts, linker)

    nt_m_indexes = get_indexes(aa, m_name, m_part, 3)
    aa_m_indexes = get_indexes(aa, m_name, m_part, 1)

    m_indexes = aa_m_indexes

    legend = "leader seqeunce | Linker sequence | cdr3 sequence | v region | j region | c region | Start | End"
    l_name = []
    l_part = []
    for i in legend.split(' | '):
        l_name.append('_'+i)
        l_part.append(i)
    l_indexes = get_indexes(legend, l_name, l_part)

    sg.theme('DarkBlue3')
    font1 = ('Courier New', 10)
    font2 = ('Courier New', 10, 'bold')
    sg.set_options(font=font1)

    layout = [
        [sg.Multiline(aa, size=(100, 10), key='-Multiline', disabled=True)],
        [sg.Multiline(legend, size=(100, 1), key='-Legend', disabled=True)],
        [sg.Push(), sg.Button('Highlight'), sg.Button('Exit'), sg.Button('NT'), sg.Button('AA', disabled=True)],
    ]

    window = sg.Window('Sequence Display', layout, finalize=True)
    multiline = window['-Multiline']
    legend = window['-Legend']
    m_widget = multiline.Widget
    l_widget = legend.Widget

    mode = 1
    while True:
        event, values = window.read()
        if event in (sg.WIN_CLOSED, 'Exit'):
            break
        elif event == 'Highlight':
            m_widget = get_highlights(m_widget, m_indexes, font2)
            l_widget = get_highlights(l_widget, l_indexes, font2)
            window['-Multiline'].update(disabled=True)
            window['-Legend'].update(disabled=True)
            window['Highlight'].update(disabled=True)
        elif event == 'NT':
             window['-Multiline'].update(nt)
             window['NT'].update(disabled=True)
             window['AA'].update(disabled=False)
             window['Highlight'].update(disabled=False)
             m_indexes = nt_m_indexes
        elif event == 'AA':
             window['-Multiline'].update(aa)
             window['NT'].update(disabled=False)
             window['AA'].update(disabled=True)
             window['Highlight'].update(disabled=False)
             m_indexes = aa_m_indexes

    window.close()


def main():
    print("Please use the appropriate 'stitchr', 'thimble', 'gui_stitchr' or 'stitchrdl' command.")
