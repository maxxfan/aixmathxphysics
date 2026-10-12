# AI × Math × Physics

Notes and experiments at the intersection of AI, mathematics and physics.

Live site: https://maxxfan.github.io/aixmathxphysics/

## Layout

- `index.html` — home page and list of papers
- `papers/<slug>/index.html` — web edition of a paper (math typeset client-side with KaTeX), with the PDF beside it
- `assets/` — shared stylesheet and the KaTeX setup (`math.js` holds LaTeX macros such as `\code`)

## Development

The site is plain static HTML served by GitHub Pages from the `main` branch root; there's no build step.
Preview it locally with:

```bash
python3 -m http.server 8000
```

Then open http://localhost:8000. Pushing to `main` redeploys the site automatically.
