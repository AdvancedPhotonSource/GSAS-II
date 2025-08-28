This is a branch for working on logging & replay of actions in GSAS-II

Changes made so far are in two sections. 

G2EventLogger: Captures wx events and logs them in such a way that they can be reproduced in a separate process (widgets are identified by strings and numbers, not by Id's or objects. 

Problem: how to track events from modal windows? It does not appear that OS-supplied ones create events and even for locall-written ones, more needs to be done to capture the input. Does the event logger need to be rgistered for all dialogs?

Also, some widgets (checkbuttons) do not appear to generate events. 

The eventlogger is currently invoked at startup, but needs to be invoked by a menu command. 

A button is provided to invoke test commands to try out different aspects of "playback."

I envision two ways this could work. GUI-scripting, where the goal is to test that the GUI is working as expected and "Demo" mode, where animation is used to show users how to do something in the GUI. Here I try out different things that could be of use.  

Note that screen coordinates need to be generated relative to the window location -- except for the menu on the Mac. 
