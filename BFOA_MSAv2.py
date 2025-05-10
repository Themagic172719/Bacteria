from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy as np
import time

class BacteriaSimulation:
    def __init__(self, path, inicio, num_bacteria=5, num_random_bacteria=1, iterations=30, gaps=1):
        self.poblacion = [bacteria(path) for _ in range(num_bacteria)]
        self.path = path
        self.inicio = inicio
        self.num_random_bacteria = num_random_bacteria
        self.iterations = iterations
        self.gaps = gaps
        self.chemio = chemiotaxis()
        self.veryBest = bacteria(path)
        self.original = bacteria(path)
        self.global_nfe = 0

        self.d_attr = 0.8
        self.w_attr = 0.5
        self.h_rep = self.d_attr
        self.w_rep = 12

    def clona_best(self, best):
        """Clone the best bacteria to the veryBest."""
        self.veryBest.matrix.seqs = np.array(best.matrix.seqs)
        self.veryBest.blosumScore = best.blosumScore
        self.veryBest.fitness = best.fitness
        self.veryBest.interaction = best.interaction

    def valida_secuencias(self):
        """Validate sequences against the original bacteria."""
        temp_bacteria = bacteria(self.path)
        temp_bacteria.matrix.seqs = np.array(self.veryBest.matrix.seqs)

        for i in range(len(temp_bacteria.matrix.seqs)):
            temp_bacteria.matrix.seqs[i] = temp_bacteria.matrix.seqs[i].replace("-", "")
            if temp_bacteria.matrix.seqs[i] != self.original.matrix.seqs[i]:
                print("*****************Secuencias no coinciden********************")
                return

    def run_simulation(self):
        """Run the bacteria simulation for the specified number of iterations."""
        for _ in range(self.iterations):
            for bacteria in self.poblacion:
                bacteria.tumboNado(self.gaps)
                bacteria.autoEvalua()

            self.chemio.doChemioTaxis(self.poblacion, self.d_attr, self.w_attr, self.h_rep, self.w_rep)
            self.global_nfe += self.chemio.parcialNFE
            
            best = max(self.poblacion, key=lambda x: x.fitness)
            if (self.veryBest is None) or (best.fitness > self.veryBest.fitness):
                self.clona_best(best)
            if _ ==self.iterations-1:
                nombre_archivo = "resultados.csv"
                fin = time.time()
                self.agregar_linea_al_archivo(nombre_archivo, f"{self.veryBest.fitness},{self.global_nfe},{fin-self.inicio}")
                print(f"Interacción: {self.veryBest.interaction}, Fitness: {self.veryBest.fitness}, NFE: {self.global_nfe}")
            
            self.chemio.eliminarClonar(self.path, self.poblacion)
            self.chemio.insertRamdomBacterias(self.path, self.num_random_bacteria, self.poblacion)
            print(f"Población: {len(self.poblacion)}")

        self.veryBest.showGenome()
        self.valida_secuencias()
    def agregar_linea_al_archivo(self, nombre_archivo, texto):
        with open(nombre_archivo, 'a') as archivo:
            archivo.write(texto + '\n')

if __name__ == "__main__":
    inicio = time.time()
    simulation = BacteriaSimulation("multiFasta.fasta", inicio)
    simulation.run_simulation()
