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
