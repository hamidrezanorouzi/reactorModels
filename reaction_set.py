from reaction import reaction

class reaction_set:
    
    def __init__(self, reactions):
        self.reactions = reactions
        
    def append_reaction(self, react : reaction):
        self.reactions.append(react)
        
    def num_reactions(self):
        return len(self.reactions)
    
    def all_net_productions(self, C, T):         # kmol/m^3.sec
       net_productions_sum = sum([reaction.net_production(C, T) for reaction in self.reactions])
       return net_productions_sum
        
    def Hrxn_dot_rate(self, T, C):               # j/m^3.sec
        return [reaction.Hrxn(T) * reaction.rate_of_reaction(C, T) for reaction in self.reactions]
