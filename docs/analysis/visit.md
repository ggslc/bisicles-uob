# Viewing BISICLES output with VisIt

VisIt is a powerful 3D visualization tool for scientific data. It can be downloaded from [VisIt's official website](https://wci.llnl.gov/codes/visit/).

## Overview

This document describes how to visualize BISICLES output using VisIt, with a focus on creating a video of ice sheet retreat over time. The example uses data from the [Pine Island Glacier example](../examples/pineisland.md).

**Before starting**: Ensure you have run the Pine Island Glacier example and have VisIt installed on a workstation with access to the plot files. You can also use VisIt's remote visualization features to view files on a remote machine.

## Opening a sequence of files

### Step 1: Browse and open the plot files

In VisIt's main window, select **File → Open**. Use the file dialog to navigate to the directory containing your plot files.

**Tip**: Use the filter `plot*hdf5` to display only plot files and hide checkpoint files (`chk*hdf5`) that VisIt cannot read.

**Smart File Grouping**: Select "Smart File Grouping" to automatically group related plot files into a single object (e.g., `plot.pigv5.1km.l1l2.4lev.*hdf5 database`).

Select this grouped object—it should appear as the "Active Source" in VisIt's main window.

### Step 2: Add visualization plots

Create the visualization by adding the following plots (accessible from the **Plots → Add** menu):

#### 1. Basal friction (pseudocolor)

1. **Plots → Add → Pseudocolor → basal_friction**
2. Double-click the "Pseudocolor - basal friction" item in the Plots list
3. In the dialog that opens:
   - Set Maximum value: `1`
   - Change Color table to: `gray`
   - Set Opacity to: `25%`

This creates a translucent overlay showing basal friction values.

#### 2. Velocity magnitude (pseudocolor, log scale)

1. **Plots → Add → Pseudocolor → Vel_magnitude**
2. Double-click the "Vel_magnitude" item in the Plots list
3. In the dialog that opens:
   - Enable log scale
   - Set Minimum value: `1`
   - Set Maximum value: `1e3` (1000 m/a)
   - Change Color table to: `bluehot`

This displays ice flow speed on a logarithmic scale.

#### 3. Mesh levels (subset with wireframe)

1. **Plots → Add → Subset → levels**
2. Double-click the "Subset - levels" item in the Plots list
3. Check the "Wireframe" box in the resulting dialog

This shows the adaptive mesh refinement (AMR) grid as colored boxes, with finer resolution represented by smaller boxes.

### Step 3: Render the visualization

Click the **Draw** button. You should now see a visualization showing:
- **Ice flow speed** (color map from the velocity magnitude plot)
- **Ice shelf outline** (translucent gray region from basal friction)
- **Mesh structure** (progressively smaller colored boxes showing AMR refinement levels, with the finest boxes at 250m resolution)

If you only see large boxes and not the refinement structure:
- Ensure "Apply subset selections to all plots" is checked (bottom of the main window)
- Use **Controls → Subset** to select all refinement levels

## Animating time-dependent data

To create an animation showing ice sheet evolution over time:

- Use the **Time slider** in the main window to step through time steps
- Click the **Playback button** (below the time slider) to animate
- Or use the playback button in the image window ("Window 1") to play the animation

As you advance through time, you'll see:
- The grounding line retreat up Pine Island Glacier's trunk
- The mesh automatically refining/coarsening to maintain high resolution at the grounding line
- Changes in velocity and basal friction patterns

## Optional enhancements

- **Add velocity vectors**: From the **Plots → Add** menu, select **Vectors** to overlay ice velocity vectors on your visualization
- **Adjust color ranges**: Experiment with different color tables and value ranges to highlight features of interest
- **Save snapshots/movies**: Use **File → Save/Export** to save individual frames or create a movie file
