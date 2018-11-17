import error as e
import application_galerkin_potential as ap

xs = [x for x, y in ap.pde_data]
ys = [y for x, y in ap.pde_data]


def analytical_function(y):
    return 18*y-8


analytical = [analytical_function(y) for x, y in ap.pde_data]
print('Valor de y: ', list(y for [x, y] in ap.pde_data))

collocation_error = e.absolute_error(analytical, ap.answer)
print('Analitico:', analytical)
print('Aproximado: ', ap.answer)
print('Erro: ', collocation_error)
print(max(collocation_error))