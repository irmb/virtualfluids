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

    return (1 / settings.viscosity
            * (- pressure_grad + settings.density * settings.force)
            * z / 2 * (settings.height - z))


def poiseuille_at_heights(settings: PoiseuilleSettings, heights):
    return [poiseuille_at_z(settings, z) for z in heights]
