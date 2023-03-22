Running hst1pass with artificial star planting using xym2mat and xym2bar
========================================================================

Instructions on how to run artificial star planting using hst1pass, combined with xym2mat and xym2bar (in lieu of hst2collate).

Software dependencies
---------------------
This runs as of March 2023 using:
     - Jay’s hst1pass code (July 2022)
     - Jay’s xym2mat and xym2bar codes (March 2020)
     - hst_flystar scripts (March 2023 version)

Conceptual overview
-------------------
First, you reduce the data as usual, using hst1pass to extract positions, then xym2mat and xym2bar to iteratively collate and average the lists to get positions in the best possible reference frame.
Next, you add artificial stars and try to recover them with hst1pass. To add the artificial stars in a consistent manner, you have to specify the frame in which you add them.
This can be done by using the MAT file of the final transformation derived from the real stars.
After extracting the positions, simply collate the output lists together with xym2bar, using the IN.xym2bat file of the final averaging derived from the real stars. You do not have to iteratively run xym2mat or xym2bar, because you are already in the correct reference frame.
All that has to be done is collate the lists.
Note also: hst1pass artificial stars planting works differently than KS2... all artificial stars are recovered SIMULTANEOUSLY, not one-by-one.
So you cannot plant stars in the same position, or too near each other.

Specific commands to run the example
------------------------------------
Run the commands in this order and in the directories.
It assumes you have the 4 flt images in 00.DATA and the fake starlist FAKE.UVW in 02.ARTSTAR.

In 01.XYM:
   - bash go1_hst1pass
   - python go2_xym2matbar.py
     
In 02.ARTSTAR:
   - bash go1_copy_mat
   - bash go2_hst1pass
   - bash go3_bar
