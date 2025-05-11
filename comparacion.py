import matplotlib.pyplot as plt
import numpy as np

# Configuración
iteraciones = np.arange(1, 31)
plt.style.use('ggplot')

# Datos (ajusta según tus resultados reales)
datos = {
    'original_fitness': np.linspace(7000, 8500, 30) + np.random.normal(0, 150, 30),
    'mejorado_fitness': np.linspace(7500, 9200, 30) + np.random.normal(0, 100, 30),
    'original_time': np.linspace(10, 150, 30),
    'mejorado_time': np.linspace(8, 100, 30)
}

# Gráfico 1: Fitness
plt.figure(figsize=(10,5))
plt.plot(iteraciones, datos['original_fitness'], 'b-', label='Original')
plt.plot(iteraciones, datos['mejorado_fitness'], 'r-', label='Mejorado')
plt.title('Comparación de Fitness')
plt.xlabel('Iteraciones')
plt.ylabel('Puntuación')
plt.legend()
plt.show()