from dataclasses import dataclass


@dataclass
class PoiseuilleSettings:
    density = 1
    viscosity = 0.005
    height = 10
    length = 1
    pressure_in = 0
    pressure_out = 0
    force = 0


def poiseuille_at_z(settings: PoiseuilleSettings, z: float):
    pressure_grad = ((settings.pressure_out - settings.pressure_in) / settings.length)

    return ((1 / settings.viscosity)
            * (- pressure_grad + settings.density * settings.force)
            * (z / 2)
            * (settings.height - z))


def poiseuille_at_heights(settings: PoiseuilleSettings, heights):
    return [poiseuille_at_z(settings, z) for z in heights]


def reynolds_number(settings: PoiseuilleSettings):
    max_v = poiseuille_at_z(settings, settings.height / 2)
    return max_v * settings.height / settings.viscosity


if __name__ == '__main__':
    sim_settings = PoiseuilleSettings()

    sim_settings.force = 2e-7
    sim_settings.viscosity = 1e-3
    sim_settings.height = 16
    print(f"v_max = ", poiseuille_at_z(sim_settings, sim_settings.height / 2))
    print(f"Re =", reynolds_number(sim_settings))

    sim_settings.viscosity *= 2
    sim_settings.height *= 2
    sim_settings.force /= 2
    print(f"v_max = ", poiseuille_at_z(sim_settings, sim_settings.height / 2))
    print(f"Re =", reynolds_number(sim_settings))

    sim_settings.viscosity *= 2
    sim_settings.height *= 2
    sim_settings.force /= 2
    print(f"v_max = ", poiseuille_at_z(sim_settings, sim_settings.height / 2))
    print(f"Re =", reynolds_number(sim_settings))
