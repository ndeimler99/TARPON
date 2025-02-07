#!/usr/bin/env python3


import cairo
import argparse

def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--telo_sequences", required=True)
    parser.add_argument("--repeat", required=True)
    parser.add_argument("--mutant", required=True)
    parser.add_argument("--telo_plot", required=True)

    return parser

def hamming_distance(seq, repeat):
    hamming = 0
    for i in range(0, len(seq)):
        if seq[i] != repeat[i]:
            hamming += 1
    return hamming

def visualize_telo_seq(telo_sequences, repeat, mutant, file_path, image_width=1500, image_height=1000):

    max_len = 0
    for seq in telo_sequences:
        if len(seq) > max_len:
            max_len = len(seq)

    telo_sequences = sorted(telo_sequences, key=len, reverse=True)

    print(max_len)
    # create canvas
    surface = cairo.ImageSurface(cairo.FORMAT_RGB24, image_width, image_height)
    ctx = cairo.Context(surface)
    ctx.rectangle(0, 0, image_width, image_height)
    ctx.set_source_rgb(1,1,1)
    ctx.fill()

    x_offset_left = 100
    x_offset_right = 50
    y_offset_top = 50
    y_offset_bottom = 50

    ctx.set_font_size(20)
    ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

    # create plot lines
    ctx.set_source_rgb(0,0,0)
    ctx.rectangle(x_offset_left - 3, y_offset_top, 3, image_height-y_offset_bottom-y_offset_top + 3)
    ctx.rectangle(x_offset_left, image_height-y_offset_bottom, image_width-x_offset_left - x_offset_right, 3)
    ctx.fill()

    # establish individual read height and nucl width
    sequence_height = (image_height-y_offset_bottom-y_offset_top) / len(telo_sequences)
    nucl_width = (image_width - x_offset_left - x_offset_right) / max_len


    ctx.set_source_rgb(0.11764705882352941,0.5333333333333333,0.8980392156862745)
    ctx.rectangle(x_offset_left + 20, y_offset_top - 30, 15, 15)
    ctx.move_to(x_offset_left + 40, y_offset_top - 16)
    ctx.text_path(repeat)
    ctx.fill()

    ctx.set_source_rgb(1,0.7568627450980392,0.027450980392156862)
    ctx.rectangle(x_offset_left + 320, y_offset_top - 30, 15, 15)
    ctx.move_to(x_offset_left+340, y_offset_top - 16)
    ctx.text_path("One Nucleotide Substitutions of {}".format(repeat))
    ctx.fill()

    if mutant != "false":
        ctx.set_source_rgb(0.8470588,0.10588235294117647,0.3764705882352941)
        ctx.rectangle(x_offset_left + 720, y_offset_top - 30, 15, 15)
        ctx.move_to(x_offset_left + 740, y_offset_top - 16)
        ctx.text_path(mutant)
        ctx.fill()
    # draw tick marks and label sizes

    ctx.set_source_rgb(0,0,0)
    for i in range(0, max_len-1000, 1000):
        if i % 2000 == 0:
            a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
            ctx.move_to(x_offset_left + i*nucl_width - width/2 + 1000*nucl_width - 2.5, image_height - y_offset_bottom + 25)
            ctx.text_path('{}'.format(i))
            ctx.fill()
        ctx.rectangle(x_offset_left + i*nucl_width - 2.5, image_height-y_offset_bottom, 5, 5)
        ctx.fill()

    for i in range(0, len(telo_sequences), 100):
        if i % 200 == 0:
            a,b,width,height,c,d = ctx.text_extents('{}'.format(i))
            ctx.move_to(x_offset_left - width - 10, image_height - y_offset_bottom - i*sequence_height + height/2)
            ctx.text_path('{}'.format(i))
        ctx.rectangle(x_offset_left-5, image_height - y_offset_bottom - i*sequence_height - 2.5, 5, 5)
        ctx.fill()

    x_rect = 100
    y_rect = y_offset_top
    height=200
    for seq in telo_sequences:
        i = 0
        while i + len(repeat) <= len(seq):
            if seq[i:i+len(repeat)] == repeat:
                i += len(repeat)
                ctx.set_source_rgb(0.11764705882352941,0.5333333333333333,0.8980392156862745)
                ctx.rectangle(x_rect, y_rect, len(repeat)*nucl_width, sequence_height)
                ctx.fill()
                x_rect += len(repeat)*nucl_width  
            elif mutant != "false" and seq[i:i+len(mutant)] == mutant:
                    i += len(mutant)
                    ctx.set_source_rgb(0.8470588,0.10588235294117647,0.3764705882352941)
                    ctx.rectangle(x_rect, y_rect, len(mutant)*nucl_width, sequence_height)
                    ctx.fill()
                    x_rect += len(mutant) * nucl_width
            elif hamming_distance(seq[i:i+len(repeat)], repeat) == 1:
                i += len(repeat)
                ctx.set_source_rgb(1,0.7568627450980392,0.027450980392156862)
                ctx.rectangle(x_rect, y_rect, len(repeat)*nucl_width, sequence_height)
                ctx.fill()
                x_rect += len(repeat)*nucl_width  
            else:
                x_rect += 1*nucl_width
                i += 1
        y_rect += sequence_height
        x_rect = 100

    surface.write_to_png(file_path)

def main(args):
    telo_seqs = []
    with open(args.telo_sequences, "r") as telo_fh:
        for line in telo_fh:
            telo_seqs.append(line.strip())

    visualize_telo_seq(telo_seqs, args.repeat, args.mutant, args.telo_plot)

if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)

