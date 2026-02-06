import numpy as np
import pytest
import sympy

from lr_reduction.user_defined_function import UserDefinedFunction, UserDefinedFunctionError, parse_user_expression


class TestParseUserExpression:
    """Tests for parse_user_expression function."""

    def test_parse_simple_formula(self):
        """Test parsing a simple formula with no parameters."""
        formula, params = parse_user_expression("formula=x+1")
        assert formula == "x+1"
        assert params == {}

    def test_parse_formula_with_parameters(self):
        """Test parsing formula with multiple parameters."""
        formula, params = parse_user_expression("formula=a*x+b, a=2.5, b=1.0")
        assert formula == "a*x+b"
        assert params == {"a": 2.5, "b": 1.0}

    def test_parse_with_name_field(self):
        """Test that name field is skipped."""
        formula, params = parse_user_expression("name=MyFunc, formula=x**2, c=3.0")
        assert formula == "x**2"
        assert params == {"c": 3.0}

    def test_parse_with_whitespace(self):
        """Test parsing with extra whitespace."""
        formula, params = parse_user_expression("  formula = x*2  ,  param = 5.5  ")
        assert formula == "x*2"
        assert params == {"param": 5.5}

    def test_parse_missing_formula_raises_error(self):
        """Test that missing formula raises ValueError."""
        with pytest.raises(ValueError, match="Missing Formula"):
            parse_user_expression("a=1.0, b=2.0")

    def test_parse_empty_parameters(self):
        """Test parsing formula with no parameters."""
        formula, params = parse_user_expression("formula=sin(x)")
        assert formula == "sin(x)"
        assert params == {}


class TestUserDefinedFunction:
    """Tests for UserDefinedFunction class."""

    def test_simple_linear_function(self):
        """Test evaluation of a simple linear function."""
        udf = UserDefinedFunction("formula=2*x+1")
        result = udf(3)
        assert result == 7

    def test_function_with_parameters(self):
        """Test function with parameters."""
        udf = UserDefinedFunction("formula=a*x+b, a=2.0, b=3.0")
        result = udf(5)
        np.testing.assert_almost_equal(result, 13.0)

    def test_function_with_numpy_array(self):
        """Test function evaluation with numpy array input."""
        udf = UserDefinedFunction("formula=x**2")
        x = np.array([1, 2, 3, 4])
        result = udf(x)
        np.testing.assert_array_equal(result, np.array([1, 4, 9, 16]))

    def test_function_with_multiple_parameters(self):
        """Test function with multiple parameters."""
        udf = UserDefinedFunction("formula=a*x**2+b*x+c, a=1.0, b=2.0, c=3.0")
        result = udf(2)
        np.testing.assert_almost_equal(result, 11.0)

    def test_function_with_trigonometric_expression(self):
        """Test function with trigonometric expressions."""
        udf = UserDefinedFunction("formula=sin(x)")
        result = udf(0)
        np.testing.assert_almost_equal(result, 0.0)

    def test_function_with_exponential(self):
        """Test function with exponential expression."""
        udf = UserDefinedFunction("formula=exp(scale*x), scale=1.0")
        result = udf(0)
        np.testing.assert_almost_equal(result, 1.0)

    def test_invalid_formula_raises_error(self):
        """Test that invalid formula raises an error."""
        with pytest.raises(sympy.SympifyError):
            UserDefinedFunction("formula=x+++")

    def test_parameter_name_mismatch_raises_error(self):
        """Test that parameter names in formula must match those provided."""
        with pytest.raises(UserDefinedFunctionError, match="Parameter mismatch"):
            UserDefinedFunction("formula=a*x+b, a=2.0, c=3.0")
