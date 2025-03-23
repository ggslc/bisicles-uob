
# Purpose

Demonstrates user of amr.evolve_ice_frac2 = true

Ice area fraction $f$ evolves according to the non-conservative advection equation
$$ \frac{\partial f}{\partial t} + (\vec{u} - \vec{u}_c) . \nabla f $$

- $f$ - fractional ice coverered area within cell
- $\vec{u}$: ice velocity
- $\vec{u}_c$: calving vector.

# Simulations

- Straight channel, 128 km $\times$ 16 km domain
- No-slip conditions at $x = 0$, $y = 0$, $y = 16 km$
- Initial state provides a calving front at $x = 96 km$
- Front held steady for $ 0 < t \leq 8 $ years by setting calving rate $u_c = u$
- Front retreats when $ 8 < t < 256 $, three cases:
	- rn: $\vec{u}_c = \vec{u} + \alpha \nabla f$ - 
	- ri: $\vec{u}_c = \vec{u} + \alpha \vec{u} / |\vec{u}|$ 
        - rp: $\vec{u}_c = \gamma \vec{u}$
- Each case repeated over 0-4 levels of refinement

# How to run

1. mk_inputs.sh \# produces inputs.\* files for each case
2. mk_ouputs.sh \# run BISICLES for each case. runs all jobs in the background, as simultaneous serial jobs
3. python3 plot_channel_frac.py - \# produces the plots. Needs libamrfile to read the hdf5 files

#Results

## Ice area and calving zone area


![image](retreat_area_time.png)

##Snaphot images

### Set retreat rate, directed along $u$

$u_c = \alpha u/|u|$ 

![image](snapshots_ri_AMR0.png)
![image](snapshots_ri_AMR1.png)
![image](snapshots_ri_AMR2.png)
![image](snapshots_ri_AMR3.png)
![image](snapshots_ri_AMR4.png)


### Set retreat rate, directed along $\nabla f$

$u_c = \alpha \nabla f$

![image](snapshots_rn_AMR0.png)
![image](snapshots_rn_AMR1.png)
![image](snapshots_rn_AMR2.png)
![image](snapshots_rn_AMR3.png)
![image](snapshots_rn_AMR4.png)

### Retreat rate proportional to $u$

$u_c = \alpha u$

![image](snapshots_rp_AMR0.png)
![image](snapshots_rp_AMR1.png)
![image](snapshots_rp_AMR2.png)
![image](snapshots_rp_AMR3.png)
![image](snapshots_rp_AMR4.png)

