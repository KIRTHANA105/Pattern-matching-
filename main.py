import streamlit as st

def naive_search(txt, pat):
    n = len(txt)
    m = len(pat)
    indices = []
    for i in range(n - m + 1):
        match = True
        for j in range(m):
            if txt[i + j] != pat[j]:
                match = False
                break
        if match:
            indices.append(i)
    return indices, "O((n - m + 1) * m)", "O(1)"

def rk_search(txt, pat, d=256, q=101):
    n = len(txt)
    m = len(pat)
    h = 1
    p = 0  # hash for pattern
    t = 0  # hash for text window
    indices = []

    for i in range(m - 1):
        h = (h * d) % q

    for i in range(m):
        p = (d * p + ord(pat[i])) % q
        t = (d * t + ord(txt[i])) % q

    for i in range(n - m + 1):
        if p == t:
            if txt[i:i + m] == pat:
                indices.append(i)
        if i < n - m:
            t = (d * (t - ord(txt[i]) * h) + ord(txt[i + m])) % q
            if t < 0:
                t += q

    return indices, "Average: O(n + m), Worst: O(n * m)", "O(1)"

def compute_lps(pat):
    m = len(pat)
    lps = [0] * m
    length = 0
    i = 1
    while i < m:
        if pat[i] == pat[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1
    return lps

def kmp_search(txt, pat):
    n = len(txt)
    m = len(pat)
    lps = compute_lps(pat)
    i = j = 0
    indices = []

    while i < n:
        if pat[j] == txt[i]:
            i += 1
            j += 1
        if j == m:
            indices.append(i - j)
            j = lps[j - 1]
        elif i < n and pat[j] != txt[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1
    return indices, "O(n + m)", "O(m)"

# Streamlit UI
st.title("Pattern Matching: Naive, Rabin-Karp, and KMP Algorithms")
st.write("This tool demonstrates how different pattern matching algorithms work and compares their results and complexities.")

txt = st.text_area("Text Input", "ABABDABACDABABCABAB")
pat = st.text_input("Pattern to Search", "ABABCABAB")

if st.button("Search"):
    if not txt or not pat:
        st.warning("Please enter both the text and the pattern.")
    elif len(txt) <= len(pat):
        st.warning("Text length should be greater than pattern length.")
    else:
        st.header("Results")

        # Naive Algorithm
        res_naive, tc_naive, sc_naive = naive_search(txt, pat)
        st.subheader("Naive Algorithm")
        st.write("Match Indices:", res_naive or "No matches found")
        st.write("Time Complexity:", tc_naive)
        st.write("Space Complexity:", sc_naive)

        # Rabin-Karp Algorithm
        res_rk, tc_rk, sc_rk = rk_search(txt, pat)
        st.subheader("Rabin-Karp Algorithm")
        st.write("Match Indices:", res_rk or "No matches found")
        st.write("Time Complexity:", tc_rk)
        st.write("Space Complexity:", sc_rk)

        # KMP Algorithm
        res_kmp, tc_kmp, sc_kmp = kmp_search(txt, pat)
        st.subheader("Knuth-Morris-Pratt (KMP) Algorithm")
        st.write("Match Indices:", res_kmp or "No matches found")
        st.write("Time Complexity:", tc_kmp)
        st.write("Space Complexity:", sc_kmp)

        st.markdown("---")
        st.caption("All algorithms return the starting indices where the pattern is found in the given text.")
