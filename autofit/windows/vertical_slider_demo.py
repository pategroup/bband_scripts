import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import AxesWidget, Slider, Button
from matplotlib.transforms import Affine2D


class VSlider(AxesWidget):
    """
    A slider representing a floating point range

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *hline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None, slidermax=None,
                 dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*

        *valinit*
            The slider initial position

        *label*
            The slider label

        *valfmt*
            Used to format the slider value

        *closedmin* and *closedmax*
            Indicate whether the slider interval is closed

        *slidermin* and *slidermax*
            Used to constrain the value of this slider to the values
            of other sliders.

        additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...)
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin,valinit,0,1, **kwargs)

        self.hline = ax.axhline(valinit,0,1, color='r', lw=1)


        self.valfmt=valfmt
        ax.set_yticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_xticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.02, label, transform=ax.transAxes,
                             verticalalignment='bottom',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.02, valfmt%valinit,
                               transform=ax.transAxes,
                               verticalalignment='top',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active  = False

    def _update(self, event):
        'update the slider position'
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event')
             or (event.name == 'button_press_event' and event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt%val)
        if self.drawon: self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson: return
        for cid, func in self.observers.iteritems():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        'remove the observer with connection id *cid*'
        try: del self.observers[cid]
        except KeyError: pass

    def reset(self):
        "reset the slider to the initial value if needed"
        if (self.val != self.valinit):
            self.set_val(self.valinit)



ax = plt.subplot(111)
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 1.0, 0.001)
a0 = 5
f0 = 3
s = a0*np.sin(2*np.pi*f0*t)
l, = plt.plot(t,s, lw=2, color='red')
plt.axis([0, 1, -10, 10])

axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axamp  = plt.axes([0.10, 0.2, 0.03, 0.65], axisbg=axcolor)

sfreq = Slider(axfreq, 'Freq', 0.1, 30.0, valinit=f0)
samp = VSlider(axamp, 'Amp', 0.1, 10.0, valinit=a0)

#axamp.set_transform(tr)
#samp.poly.set_transform(tr)

def update(val):
    amp = samp.val
    freq = sfreq.val
    l.set_ydata(amp*np.sin(2*np.pi*freq*t))
    plt.draw()
sfreq.on_changed(update)
samp.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    sfreq.reset()
    samp.reset()
button.on_clicked(reset)


plt.show()

