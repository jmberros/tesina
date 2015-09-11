from IPython.display import display, Math


def num_format(x):
    return '{:,.0f}'.format(x)


def format_numbers(df):
    """Returns a copy of the DataFrame with its numbers formatted"""
    df = df.copy()

    for series in df:
        if df[series].dtype in ['float64', 'int64']:
            df[series] = df[series].apply(num_format)

    return df


def m(math_text):
    """Returns a math LaTeX text fot IPython"""
    return display(Math(math_text))
