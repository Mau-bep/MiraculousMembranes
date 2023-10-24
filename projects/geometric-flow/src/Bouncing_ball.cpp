#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <iostream>
#include <string>
#include <unistd.h> // for usleep function

int main()
{
    // open display and create window
    Display* display = XOpenDisplay(NULL);
    int screen = DefaultScreen(display);
    Window window = XCreateSimpleWindow(display, RootWindow(display, screen),
        0, 0, 400, 400, 0, BlackPixel(display, screen), WhitePixel(display, screen));

    // set window properties
    XSelectInput(display, window, ExposureMask);
    XMapWindow(display, window);

    // create graphics context
    GC gc = XCreateGC(display, window, 0, NULL);

    // set initial ball position and velocity
    int ball_x = 100;
    int ball_y = 150;
    int vel_x = -2;
    int vel_y = 3;

    // set initial ball color
    XSetForeground(display, gc, WhitePixel(display, screen));

    // loop until 'q' key is pressed
    while (true)
    {
        // clear window
        XClearWindow(display, window);

        // check for collisions with edges
        if (ball_x <= 10 || ball_x >= 370)
        {
            vel_x = -vel_x;

            // change ball color
            std::string color = (ball_y % 2 == 0) ? "red" : "green";
            XColor xcolor;
            Colormap cmap = DefaultColormap(display, screen);
            XParseColor(display, cmap, color.c_str(), &xcolor);
            XAllocColor(display, cmap, &xcolor);
            XSetForeground(display, gc, xcolor.pixel);
        }
        if (ball_y <= 10 || ball_y >= 350)
        {
            vel_y = -vel_y;

            // change ball color
            std::string color = (ball_x % 2 == 0) ? "blue" : "yellow";
            XColor xcolor;
            Colormap cmap = DefaultColormap(display, screen);
            XParseColor(display, cmap, color.c_str(), &xcolor);
            XAllocColor(display, cmap, &xcolor);
            XSetForeground(display, gc, xcolor.pixel);
        }

        // move ball
        ball_x += vel_x;
        ball_y += vel_y;

        // draw ball
        XFontStruct* font_info = XLoadQueryFont(display, "fixed");
        XSetFont(display, gc, font_info->fid);
        XDrawString(display, window, gc, ball_x - 20, ball_y + 10, "DVD", 3);

        // flush changes to window
        XFlush(display);

        // pause for a short time to slow down animation
        usleep(10000);

        // check for quit key
        XEvent event;
        if (XCheckTypedEvent(display, KeyPress, &event) && event.xkey.keycode == XKeysymToKeycode(display, XK_q))
        {
            break;
        }
    }

    // clean up resources and close display
    XFreeGC(display, gc);
    XDestroyWindow(display, window);
    XCloseDisplay(display);

    return 0;
}