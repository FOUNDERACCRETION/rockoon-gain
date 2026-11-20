"""Cross-check of the trajectory integrator against RocketPy (an independent, established flight-dynamics code).

The baseline vehicle is flown from the ground and from 30 km in RocketPy with the same mass, constant vacuum thrust,
pressure-thrust correction (exit area), drag curve and inverse-square gravity, and in our own integrator.
Writes results/rocketpy_check.json."""
import json, math, os, warnings
import numpy as np
import rockoon_model as rm
from rocketpy import Environment, Rocket, Flight, GenericMotor, Function

warnings.filterwarnings('ignore')
ROOT = os.path.dirname(rm.HERE)
M0, MP, D = 301.27, 151.27, 0.40


def ours(z0):
    r = rm.vehicle_run(M0, MP, z0, d=D)
    return r


def rocketpy_run(z0, sim_mode):
    r = ours(z0)
    nz = rm.Nozzle(z0, 4 * M0 * rm.G0, D)
    Fvac = nz.mdot * rm.G0 * nz.isp_vac
    tb = MP / nz.mdot
    gfun = lambda h: rm.G0 * (rm.R_EARTH / (rm.R_EARTH + h)) ** 2
    env = Environment(gravity=gfun, elevation=z0, max_expected_height=400e3)
    # RocketPy's own standard atmosphere below 80 km (its tabulation ends near 86 km); above 80 km, where the
    # drag is negligible, the isothermal continuation of the paper's model.
    std = Environment(elevation=0.0, max_expected_height=86e3)
    std.set_atmospheric_model(type='standard_atmosphere')
    P = lambda h: float(std.pressure(h)) if h < 80e3 else rm.atmosphere(h)[1]
    T = lambda h: float(std.temperature(h)) if h < 80e3 else rm.atmosphere(h)[0]
    env.set_atmospheric_model(type='custom_atmosphere', pressure=P, temperature=T, wind_u=0.0, wind_v=0.0)
    # the Environment tabulates gravity only up to 80 km and holds it constant above; re-tabulate to 400 km
    env.max_expected_height = 400e3
    env.gravity = Function(gfun, 'height (m)', 'gravity (m/s2)').set_discrete(0.0, 400e3, 801)
    motor = GenericMotor(thrust_source=lambda t: Fvac, burn_time=tb, chamber_radius=D / 2 * 0.9, chamber_height=2.0,
                         chamber_position=1.0, propellant_initial_mass=MP, nozzle_radius=math.sqrt(nz.Ae / math.pi),
                         dry_mass=0.0, reference_pressure=0.0)
    cd = Function(lambda M: rm.cd_generic(M), inputs='Mach', outputs='Cd')
    rocket = Rocket(radius=D / 2, mass=M0 - MP, inertia=(500.0, 500.0, 5.0), power_off_drag=cd, power_on_drag=cd,
                    center_of_mass_without_motor=2.0, coordinate_system_orientation='tail_to_nose')
    rocket.add_motor(motor, position=0.0)
    if sim_mode == '6 DOF':
        rocket.add_nose(length=1.2, kind='tangent', position=4.0)
        rocket.add_trapezoidal_fins(n=4, root_chord=0.6, tip_chord=0.3, span=0.4, position=0.65, sweep_angle=30.0)
    fl = Flight(rocket=rocket, environment=env, rail_length=2.0, inclination=90.0, heading=0.0, terminate_on_apogee=True,
                max_time=1000, rtol=1e-9, atol=1e-6, simulation_mode=sim_mode)
    return dict(apogee=float(fl.apogee), v_max=float(fl.max_speed), t_apogee=float(fl.apogee_time), q_max=float(fl.max_dynamic_pressure),
                ours_apogee=float(r['apogee']), ours_vb=float(r['v_b']), ours_q=float(r['q_max']), Fvac=Fvac, tb=tb, Ae=nz.Ae)


def main():
    out = {}
    for tag, z0 in (('ground', 0.0), ('balloon', 30e3)):
        for mode in ('3 DOF', '6 DOF'):
            try:
                out[tag + '_' + mode.replace(' ', '')] = rocketpy_run(z0, mode)
            except Exception as e:  # report but do not fail the pipeline
                out[tag + '_' + mode.replace(' ', '')] = dict(error=repr(e))
    with open(os.path.join(ROOT, 'results', 'rocketpy_check.json'), 'w') as f:
        json.dump(out, f, indent=1)
    for k, v in out.items():
        print(k, v)


if __name__ == '__main__':
    main()
