import sympy as sp


class UserDefinedFunctionError(ValueError):
    pass


def parse_user_expression(s: str):
    """
    Parse function as a string into formula and params.
    If entered from a Mantid function should skip the name part.

    Parameters
    ----------
    s : str
        The user-defined function string.

    Returns
    -------
    tuple
        (formula, params):
            formula : str
                The mathematical expression as a string.
            params : dict
                A dictionary of parameter names and their float values.

    Raises
    ------
    ValueError
        If the formula is missing in the input string.
    """
    parts = [p.strip() for p in s.split(",")]

    formula = None
    params = {}

    for part in parts:
        key, value = [x.strip() for x in part.split("=", 1)]

        if key.lower() == "formula":
            formula = value
        elif key.lower() == "name":
            continue
        else:
            params[key] = float(value)

    if formula is None:
        raise ValueError("Missing Formula=...")

    return formula, params


class UserDefinedFunction:
    """Parses and evaluates a user-defined mathematical expression.

    Attributes
    ----------
    expr_str : str
        The original formula string provided by the user.

    Raises
    ------
    SympifyError
        If the formula cannot be parsed into a valid sympy expression.
    ValueError
        If there is a mismatch between parameters used in the formula and those provided.
   """
    def __init__(self, definition: str):
        self.expr_str, self.params = parse_user_expression(definition)

        self.x = sp.symbols("x")
        self.param_symbols = {k: sp.symbols(k) for k in self.params}

        self.expr = sp.sympify(
            self.expr_str,
            locals={"x": self.x, **self.param_symbols}
        )

        # Validate that all symbols in the expression are accounted for in the parameters
        expr_params = self.expr.free_symbols - {self.x}
        defined_params = set(self.param_symbols.values())

        if expr_params != defined_params:
            missing = expr_params - defined_params
            extra = defined_params - expr_params
            raise UserDefinedFunctionError(
                "Parameter mismatch in user-defined function. "
                f"Missing parameters: {[str(p) for p in missing]}. "
                f"Extra parameters: {[str(p) for p in extra]}."
            )

        self._f = sp.lambdify(
            (self.x, *self.param_symbols.values()),
            self.expr,
            modules="numpy"
        )

    def __call__(self, x):
        return self._f(x, *self.params.values())
