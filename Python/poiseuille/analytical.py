from dataclasses import dataclass


class PoiseuilleSettings:

    def __init__(self):
        self.density = 1
        self.viscosity = 0.005
        self.height = 10
        self.length = 1
        self.pressure_in = 0
        self.pressure_out = 0
        self.force = 0


def poiseuille_at_z(settings: PoiseuilleSettings, z: float):
    pressure_grad = ((settings.pressure_out - settings.pressure_in) / settings.length)

    return (1 / settings.viscosity
            * (- pressure_grad + settings.density * settings.force)
            * z / 2 * (settings.height - z))


def poiseuille_at_heights(settings: PoiseuilleSettings, heights):
    return [poiseuille_at_z(settings, z) for z in heights]


if __name__ == '__main__':
    # h1 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    # h2 = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5]
    settings = PoiseuilleSettings()
    settings.force = 1e-8
    settings.height = 32

    # print(max(poiseuille_at_heights(settings, h1)))
    # print(max(poiseuille_at_heights(settings, h2)))

    v = poiseuille_at_z(settings, 16)
    print(v)