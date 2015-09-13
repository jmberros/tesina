def annotate_bars(ax, base=1, decimals=1, suffix='', fontsize=12):
    """Annotates matplotlib bars in a barplot. Returns the texts objects."""
    texts = []

    for bar in ax.patches:
        height = bar.get_height()
        x = bar.get_x() + bar.get_width()/2  # centered text
        y = height * 1.05  # text slightly above the bar
        if decimals == 0:
            text = '{}{}'.format(int(height/base), suffix)
        else:
            text = '{}{}'.format(round(height/base, decimals), suffix)
        texts.append(ax.text(x, y, text, fontsize=fontsize, ha='center'))

    return texts


def remove_unnecessary_lists_from_df(df):
    df = df.copy()

    for field_name, s in df.iteritems():
        if not all([type(e) == list for e in s]):
            continue

        # Drop empty lists
        if all([len(e) == 0 for e in s]):
            df.drop(field_name, axis=1)

        # Remove unncessary list wrapping of [just one element]
        if all([len(e) == 1 for e in s]):
            df[field_name] = s.map(lambda e: e[0])

    return df
