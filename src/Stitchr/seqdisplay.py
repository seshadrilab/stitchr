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
    names.append("linker")
    parts.append(linker)
    names.append("Start")
    parts.append('M')
    names.append('End')
    parts.append('*')
    return zip(names, parts)


def get_indexes(seq, parts):
    """
    Input: A string DNA sequence and a dictionary of gene regions and their DNA sequence
    Output: A list of tuples where each gene region occurs in DNA sequence
    """
    names = []
    index1s = []
    index2s = []
    for name, part in parts:
        if part in seq:
            if name == 'Start':
                start = seq.index(part)
                index1s.append(f'1.{start}')
                index2s.append(f'1.{start + 1}')
                names.append(name)
            else:
                for m in re.finditer(part, seq):
                    index1s.append(f'1.{m.start()}')
                    index2s.append(f'1.{m.end()}')
                    names.append(name)
    indexes = list(zip(names, index1s, index2s))
    return(indexes)


def display(seq, parts, linker):
    """
    Input: A string DNA sequence and a dictionary of gene regions and their DNA sequence
    Output: A GUI display of the DNA sequence that highlights different regions
    """
    parts = get_partlist(parts, linker)
    indexes = get_indexes(seq, parts)

    sg.theme('DarkBlue3')
    font1 = ('Courier New', 10)
    font2 = ('Courier New', 10, 'bold')
    sg.set_options(font=font1)

    layout = [
        [sg.Multiline(seq, size=(100, 20), key='-Multiline')],
        [sg.Push(), sg.Button('Highlight'), sg.Button('Exit')],
    ]

    window = sg.Window('Sequence Display', layout, finalize=True)
    multiline = window['-Multiline']
    widget = multiline.Widget
    while True:
        event, values = window.read()
        if event in (sg.WIN_CLOSED, 'Exit'):
            break
        elif event == 'Highlight':
            for name, index1, index2 in indexes:
                    # Make it so that a value is stored in in indexes that will check against defined (v, j, cdr3, c, start, stop)
                    if "_cdr3" in name:
                        widget.tag_config('BLACK', foreground='white', background='black', font=font2)
                        widget.tag_add('BLACK', index1, index2)
                    elif "_l" in name:
                        widget.tag_config('PURPLE', foreground='white', background='purple', font=font2)
                        widget.tag_add('PURPLE', index1, index2)
                    elif "linker" in name:
                        widget.tag_config('BLUE', foreground='white', background='blue', font=font2)
                        widget.tag_add('BLUE', index1, index2)
                    elif "Start" in name:
                        widget.tag_config('GREEN', foreground='white', background='green', font=font2)
                        widget.tag_add('GREEN', index1, index2)
                    elif "*" in name:
                        widget.tag_config('RED', foreground='white', background='red', font=font2)
                        widget.tag_add('RED', index1, index2)
                    elif "_v" in name:
                        widget.tag_config('ORANGE', foreground='black', background='orange', font=font2)
                        widget.tag_add('ORANGE', index1, index2)
                    elif "_c" in name:
                        widget.tag_config('PINK', foreground='white', background='brown', font=font2)
                        widget.tag_add('PINK', index1, index2)

            window['-Multiline'].update(disabled=True)
    window.close()


def main():
    print("Please use the appropriate 'stitchr', 'thimble', 'gui_stitchr' or 'stitchrdl' command.")
