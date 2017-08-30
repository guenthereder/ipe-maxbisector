# maxbisector

*maxbisector* is an *ipelet* for the drawing tool ipe. It constructs the weighted
L-infinity bisector between two points.

## Usage
* use marks as points
* use the symbol size as weight
* select two marks and find *maxbisector* in the ipelets menu

## Install
* copy or link *maxbisector* to HOME/.ipe/ipelets
   * example: `ln {maxbisector.git}/maxbisector.lua
     $HOME/.ipe/ipelets/maxbisector.lua`
* (optional) customize ipe using a
  [customize.lua](http://git.sthu.org/?p=ipestuff.git;a=blob_plain;f=customize.lua;hb=HEAD)
ipelet (linked from [Stefan Huber](https://www.sthu.org/misc/ipe.html)) and add
an additional `shortcuts.ipelet_1_maxbisector = "Ctrl+D"`

## Ipe
* [Ipe](http://ipe.otfried.org/) - extendable drawing editor

