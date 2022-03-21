from ccinput.constants import SYN_TYPES


def format_dict_str(d, name, label2="Synonyms"):
    arr_l = [len(name)]
    for k, syns in d.items():
        arr_l.append(len(k))
        arr_l += [len(s) for s in syns]

    ll = max(arr_l)

    tab = """{} ========\n{:<{}} {}\n{} ========\n""".format(
        ll * "=", name, ll, label2, ll * "="
    )

    for k, syns in sorted(d.items(), key=lambda i: i[0]):
        tab += "{:<{}} {}\n".format(k, ll, syns[0] if len(syns) > 0 else "")
        if len(syns) > 1:
            for s in syns[1:]:
                tab += "\n{} {}\n".format(ll * " ", s)

    tab += "{} ========".format(ll * "=")
    return tab.strip()


def format_dict_enum(d, name):
    return format_dict_str({k.name: v for k, v in d.items()}, name)


def format_calc_types():
    return format_dict_str(
        {j[0]: [str(i)] for i, j in SYN_TYPES.items()},
        "Calculation type",
        label2="Code",
    )
