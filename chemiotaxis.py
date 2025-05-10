import math
import random
from bacteria import bacteria

class chemiotaxis:
    def __init__(self):
        self.parcialNFE = 0  

    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 0.0
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore) ** 2
            total += d * math.exp(w * diff)
        return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
        return attract + repel

    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction

    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        self.parcialNFE = sum(bacteria.NFE for bacteria in poblacion)
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
            bacteria.NFE = 0

    def eliminarClonar(self, path, poblacion):
        """Elimina el 50% de las bacterias con menor fitness y las reemplaza con clones"""
        poblacion.sort(key=lambda x: x.fitness)
        poblacion[:] = poblacion[len(poblacion) // 2:]  # Mantener solo el 50% superior
        poblacion.extend(self.clonacion(path, poblacion))

    def clonacion(self, path, poblacion):
        best_fitness = max(bacteria.fitness for bacteria in poblacion)
        return [self._crear_clone(bacteria, path, best_fitness) for bacteria in poblacion]

    def _crear_clone(self, bacteria, path, best_fitness):
        newBacteria = bacteria.clonar(path)
        mutacion = int((best_fitness - bacteria.fitness) / 10)
        newBacteria.tumboNado(mutacion)
        newBacteria.autoEvalua()
        return newBacteria

    def randomBacteria(self, path):
        bact = bacteria(path)
        bact.tumboNado(random.randint(1, 10))
        return bact 

    def insertRamdomBacterias(self, path, num, poblacion):
        for _ in range(num):
            poblacion.append(self.randomBacteria(path))
            poblacion.sort(key=lambda x: x.fitness)
            del poblacion[0]
