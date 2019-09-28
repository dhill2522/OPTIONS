=typical pwr model -- 4 inch cold leg break 36.05 check case
*                    Configuration Control Problem
*       This problem is a simulation of a four loop presurized reactor
*  undergoing a small break.  Loop containing break is modeled as a
*  single loop but the other three loops are coalesced into one loop.
*  Modeling does not now follow all recommended modeling practices but
*  problem is still good test of many features of code.  Problem uses
*  standard matrix techniques.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*1 47 53
100 restart transnt
*
102  si   si
*   restart#  name reset
103   -1      zionpwr_PERCS_SBO_6in.r
*104  restart  zionpwr_PERCS_SBO_6in.r
105  5.0 6.0  1200.0
*******************************************************
*                                                     *
*               time step cards                       *
*                                                     *
*******************************************************
*
*     time minstep maxstep copt pfreq majed  rsrtf
201 10000.0  1.0e-5   0.1   23   10  100    100000   *1.0e-10, 1.0e-3
* 202 70000.0 1.0e-10   0.1   23   36000 36000 100000
* 203 520000.0 1.0e-10   0.1  23   1000 1000 100000
* 204 2592100.0  1.0e-10   0.1   23   40000  500000    100000
501  time   0            ge   null    0      100.0        l
*
*end
.