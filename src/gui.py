import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

class MyWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="ProteinCoLoc.jl")
        self.set_default_size(400, 200)

        # create notebook
        notebook = Gtk.Notebook()

        # create pages
        for i in range(1, 3):
            # create box for page content
            box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)

            # add some content to the box
            button = Gtk.Button(label="Load image")
            box.pack_start(button, True, True, 0)

            # create a box for tab
            tab_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)

            # create a button for tab
            button = Gtk.Button(label="Load images")
            tab_box.pack_start(button, True, True, 0)

            # create a label for tab
            label = Gtk.Label("Load images")
            tab_box.pack_start(label, True, True, 0)

            # add page to notebook
            notebook.append_page(box, tab_box)

        # add notebook to window
        self.add(notebook)

win = MyWindow()
win.connect("destroy", Gtk.main_quit)
win.show_all()
Gtk.main()