from gi.repository import Gtk
from gi.repository import Gdk
import time
from gi.repository import GObject as gobject
from math import pi, atan2, cos, sin, sqrt
#from test_verlet import WorldLawsEngine
from random import random as rand
import random
import numpy as np
import copy
from numpy.linalg import norm
import test_GEP as tgep
import _GEP as gep

SIZE = 20
WIN_SIZE = 800, 600
PI2 = pi*2

def message(ctx, pos, msg, size=30, c=(0, 0, 0)):
    ctx.set_font_size(size)
    ctx.set_source_rgb(c[0], c[1], c[2])
    ctx.move_to(pos[0], pos[1])
    ctx.show_text(msg)

def draw_centered_text(ctx, pos, msg, font_size_in_px=30, c=(0, 0, 0)):
    ctx.set_font_size(font_size_in_px * 72.0 / 96.0)
    w, h = ctx.text_extents(msg)[2:4]
    ctx.set_source_rgb(c[0], c[1], c[2])
    ctx.move_to(pos[0] - w*0.5, pos[1] + h*0.5)
    ctx.show_text(msg)

def draw_exp_tree(ctx, scr_size, layers, x_spacing, y_spacing):
    """ node = (text, children_number) """
    m = len(layers)
    for l in layers:
        if len(l) > m:
            m = len(l)

    print(scr_size)
    x_rad = (scr_size[0] - (m + 1)*x_spacing)/m
    y_rad = (scr_size[1] - (m + 1)*y_spacing)/m
    circle_rad = min(x_rad, y_rad)/2.0
    h_spacing = scr_size[1]/len(layers)
    curr_pos = [0, h_spacing*0.5]
    layers_pos = []
    text_h = circle_rad
    for l in layers:
        w_spacing = scr_size[0]/len(l)
        curr_pos[0] = w_spacing*0.5
        pos_layer = []
        for node in l:
            ctx.set_line_width(2)
            ctx.set_source_rgb(0, 0, 0)
            ctx.move_to(*curr_pos)
            ctx.new_sub_path()
            ctx.arc(curr_pos[0], curr_pos[1], circle_rad, 0, PI2)
            pos_layer.append((curr_pos[0], curr_pos[1], node[1]))
            ctx.stroke()
            draw_centered_text(ctx, curr_pos, node[0], text_h)
            curr_pos[0] += w_spacing
        layers_pos.append(pos_layer)
        curr_pos[1] += h_spacing

    for i, l in enumerate(layers_pos[:-1]):
        k = 0
        for node in l:
            for j in range(node[2]):
                ctx.set_line_width(2)
                ctx.set_source_rgb(0, 0, 0)
                ctx.move_to(node[0], node[1] + circle_rad)
                ctx.line_to(layers_pos[i + 1][k][0], layers_pos[i + 1][k][1] - circle_rad)
                k += 1
                ctx.stroke()







class Win(Gtk.Window):

    def __init__(self, size):
        super(Win, self).__init__()

        self.init_ui(size)

    def init_ui(self, size):

        self.darea = Gtk.DrawingArea()
        self.darea.connect("draw", self.redraw)
        self.add(self.darea)

        self.set_title("Much balls.")
        self.resize(size[0], size[1])
        self.size = size
        self.set_position(Gtk.WindowPosition.CENTER)
        self.connect("delete-event", Gtk.main_quit)
        self.connect('destroy', lambda w: Gtk.main_quit())
        self.connect("motion_notify_event", self.motion_notify_event)
        self.connect('key_press_event', self.key_pressed)
        #self.connect('button_press_event', self.mouse_click)
        self.show_all()
        self.pause = True
        self.last = time.time()
        gobject.timeout_add(500, self.tick) # Go call tick every 50 whatsits.

    def motion_notify_event(self, widget, event):
	    self.mx, self.my, = event.x, self.size[1] - event.y

    def key_pressed(self, event, event2):
        if event2.keyval == ord('p'):
            self.pause = not self.pause

    def redraw(self, wid, ctx):
        # paint background
        ctx.set_source_rgb(1, 1, 1) # blue
        ctx.rectangle(0, 0, self.size[0], self.size[1])
        ctx.fill()

        #message(ctx, (20, 40), 'fps: %0.1f' % (1/(time.time() - self.last)))
        self.last = time.time()
        self.draw(ctx)


    def draw(self, ctx):
        gene = "Q*+-abcd"
        gene = "Q*+*a*Qaaba"
        layers_to_draw = gep.get_code_and_args_n_from_layers(gep.parse_gene_into_tree(gene,
            tgep.get_abc()))
        draw_exp_tree(ctx, self.size, layers_to_draw, 20, 20)
        """global qqq, bbb
        ctx.save()
        ctx.scale(1, -1)
        ctx.translate(0, -self.size[1])

        ctx.restore()
        """



    #@profile
    def tick(self):
        if self.pause:
            return True
        self.darea.queue_draw()
        #IMPORTANT TO RETURN TRUE TO TIMEOUT_ADD
        #Gtk.main_quit()
        return True


def main():
    app = Win(WIN_SIZE)
    #Gtk.gtk_window_set_modal(app, True)
    Gtk.main()

if __name__ == '__main__':
    main()

