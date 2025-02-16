import numpy as np 


class ChemicalReaction:
    def __init__(self, reaction_string, reaction_rate):
        self.reaction_string = reaction_string
        self.reaction_rate = reaction_rate

        self.left_side, self.right_side = reaction_string.split("->")
        self.reactants = self.get_components(self.left_side.strip())
        self.products = self.get_components(self.right_side.strip())

    def get_components(self, components_string):
        molecules = components_string.split("+")
        molecules_map = {}

        for m in molecules:
            num = 1
            
            if m[0].isdigit():
                num = int(m[0])
                molecules_map[m[1:]] = num
            else:
                molecules_map[m] = num

        return molecules_map
    
    def __str__(self) -> str:
        return self.reaction_string
    

reactions = []
molecules_types = set() 

while True:
    reaction = input(f"Unestite {len(reactions)+1}. reakciju (ili 'exit' za izlaz): ")
    if reaction.lower() == 'exit':
        break
    
    k = input("Unesite vrednost koeficijenta reackije:")
    reactions.append(ChemicalReaction(reaction, float(k)))    

for r in reactions:
    molecules_types = molecules_types | set(r.reactants.keys())
    molecules_types = molecules_types | set(r.products.keys())


molecules_concetration = {}
molecules_types = list(molecules_types)

for m in molecules_types:
    val = float(input(f"Unesite pocetnu {m} (0 ukoliko se molekul na pocetku ne nalazi medju reaktantima):"))
    molecules_concetration[m] = val


# koncetracija molekula tokom vremena, za pocetne koncetracije uzimam da je t=0
    
def derivative(p, x_p):
    
    ret_val = np.zeros(len(x_p))
    for i, mol in enumerate(molecules_types):
        for r in reactions:
            count = 0
            reac = list(r.reactants.keys())
            prod = list(r.products.keys())

            if mol in reac:
                count -= r.reactants[mol]
            if mol in prod:
                count += r.products[mol]

            conc = 1

            for reactant in reac:
                conc *= ((x_p[molecules_types.index(reactant)])**r.reactants[reactant])

            # if len(reac) == 1:
            #     conc = x_p[molecules_types.index(reac[0])]**2
            # else:
            #     conc = x_p[molecules_types.index(reac[0])] * x_p[molecules_types.index(reac[1])] 
                
            ret_val[i] += count * conc * r.reaction_rate
    
    return ret_val
    

def calculate_solution(points_, molecules_concetration_):
    values = []
    values.append(np.array([val for val in molecules_concetration.values()]))
    h = points_[1] - points_[0]

    for i in range(len(points_)):
        if i == 0:
            continue
            
        k1 = h*derivative(points_[i-1], values[-1])
        k2 = h*derivative(points_[i-1] + h/2, values[-1] + h*k1/2)
        k3 = h*derivative(points_[i-1] + h/2, values[-1] + h*k2/2)
        k4 = h*derivative(points_[i], values[-1] + h*k3)

        values.append(values[-1] + 1/6*(k1 + 2*k2 + 2*k3 + k4))


    return values


points = np.linspace(0,1, 1000)
values = calculate_solution(points, molecules_concetration)


values = np.array(values)

import matplotlib.pyplot as plt

dim = values.shape[1]
begin_point = values[:,0]

if dim == 2:
    plt.plot(values[:, 0], values[:, 1], marker='o', linestyle='-')  # Spajanje tačaka linijom u 2D
    plt.xlabel("Concetration of " + molecules_types[0])
    plt.ylabel("Concetration of " + molecules_types[1])
  #  plt.scatter(begin_point[0], begin_point[1], color='red', marker='o', label='Begin point') 
elif dim == 3:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(values[:, 0], values[:, 1], values[:, 2], marker='o', linestyle='-')  # Spajanje tačaka linijom u 3D
   # plt.scatter(begin_point[0], begin_point[1], begin_point[2], color='red', marker='o', label='Begin point') 

    ax.set_xlabel("Concetration of " + molecules_types[0])
    ax.set_ylabel("Concetration of " + molecules_types[1])
    ax.set_zlabel("Concetration of " + molecules_types[2])


plt.show()

# iscrtavanja grafika zavisnosti koncetracija od vremena
fig, axs = plt.subplots(1, dim, figsize=(15, 5))

for i in range(dim):
    axs[i].plot(points, values[:, i], color="orange")
    axs[i].set_xlabel('Time')
    axs[i].set_ylabel(f'Concetration of {molecules_types[i]}')
    axs[i].set_title(f'Concetration of {molecules_types[i]} through time: ')

plt.tight_layout()
plt.show()

plt.plot(values[:, molecules_types.index("a")], values[:, molecules_types.index("c")], marker='o', linestyle='-')  # Spajanje tačaka linijom u 2D
plt.xlabel("Concetration of " + "a")
plt.ylabel("Concetration of " + "c")
plt.show()
  #  plt.scatter(begin_point[0], begin_point[1], color='red', marker='o', label='Begin point') 

