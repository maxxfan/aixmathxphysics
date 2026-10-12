// Typeset $…$ and $$…$$ with KaTeX. Macros mirror the LaTeX sources of our papers.
(function () {
  if (!window.renderMathInElement) return;
  window.renderMathInElement(document.body, {
    delimiters: [
      { left: "$$", right: "$$", display: true },
      { left: "$", right: "$", display: false },
    ],
    macros: {
      "\\code": "[\\![#1,#2,#3]\\!]_{#4}",
      "\\GF": "\\mathbb{F}_{#1}",
      "\\wt": "\\operatorname{wt}",
    },
    throwOnError: false,
  });
})();
