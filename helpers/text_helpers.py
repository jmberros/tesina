def trunc_text(text, max_len=25):
    if len(text) <= max_len:
        return text

    return text[:max_len] + "..."
