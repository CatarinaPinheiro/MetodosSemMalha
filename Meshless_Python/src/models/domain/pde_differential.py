class PdeDifferential:
    def __init__(self, variables):
        self.variables = variables

    @property
    def differential(self):
        if self.variables == 'x':
            return []