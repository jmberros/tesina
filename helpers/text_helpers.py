def trunc_text(text, max_len=25):
    if len(text) <= max_len:
        return text

    return text[:max_len] + "..."

def sanitize_string_for_filename(string):
    keepcharacters = (' ', '.', '_')
    sanitized = "".join(c for c in filename if c.isalnum() or c in keepcharacters)
    sanitized = sanitized.replace(" ", "_")
    return sanitized.rstrip()
