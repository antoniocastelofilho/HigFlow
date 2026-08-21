# Contributor photographs

Photographs used by the **People** section of the top-level `README.md`.

A photograph of an identifiable person is personal data. Nothing goes in this
directory without that person's explicit, informed consent — see
[Getting consent](#getting-consent) below.

## Image requirements

| | |
|---|---|
| **Aspect** | Square. Crop before committing; the README cannot crop for you |
| **Source size** | 400 × 400 px — displayed at 100 px, so this stays sharp on high-density screens |
| **Format** | `.jpg` for photographs, `.png` only if the image genuinely needs transparency |
| **File size** | **Under 100 KB.** 60–80 KB is comfortable at this resolution |
| **Filename** | `firstname-lastname.jpg`, lowercase, ASCII, hyphens for spaces — `antonio-castelo-filho.jpg` |
| **Content** | A head-and-shoulders photograph. Any image the person is happy to have on a public page is fine |

Keep the files small. The repository just shed 76 MB of tracked build output;
a dozen unoptimised photographs would quietly put a slice of it back.

To resize and compress:

```bash
# ImageMagick — crop to a centred square, resize, compress
magick input.jpg -gravity center -crop 1:1 +repage -resize 400x400 -quality 82 \
    docs/images/contributors/firstname-lastname.jpg
```

Check the result before committing:

```bash
ls -l docs/images/contributors/
```

## Adding someone

1. Obtain consent (below) and record it.
2. Add the image to this directory, meeting the requirements above.
3. Uncomment or copy one cell in the **People** section of `README.md` and fill
   in the four fields: image path, alt text, name, and role or affiliation.
4. Keep the table at **four cells per row**, so the grid stays even on a narrow
   window.

The markup for one person:

```html
<td align="center" width="25%">
  <img src="docs/images/contributors/firstname-lastname.jpg" width="100" alt="Firstname Lastname"><br>
  <sub><b>Firstname Lastname</b></sub><br>
  <sub>Role or affiliation</sub>
</td>
```

The `alt` text matters: it is what a screen reader announces and what shows if
the image fails to load. Use the person's name, not "photo" or "avatar".

## Getting consent

Ask each person directly, and be specific about what they are agreeing to. A
message that covers it:

> I am adding a section to the HigFlow README with photographs of the people who
> have contributed to the project. Would you be willing to have yours included?
>
> If so: the photograph and your name would appear on the public repository page
> at github.com/antoniocastelofilho/HigFlow, visible to anyone. You can send me a
> photograph you would like used, or tell me to take one from somewhere. You can
> ask me to remove it at any time and I will, no explanation needed.
>
> If you would rather not, that is completely fine — your name stays in the
> Authors list either way, which comes from the commit history.

Points worth keeping in that message:

- **Where it will appear and who can see it.** "Public repository page" is
  concrete; "the project" is not.
- **That they choose the photograph.** People care which picture of them is
  public.
- **That withdrawal is easy and unconditional.** This is what makes the consent
  meaningful rather than a formality.
- **That declining costs them nothing.** The Authors list is derived from git
  history and is unaffected.

Record who agreed and when, somewhere outside this repository — an email thread
or a message log is enough. If someone asks to be removed, delete the image and
the table cell in the same commit.

## Who is in the commit history

From `git shortlog -sne --all`, with `.mailmap` applied. This is the set of
people to ask; it is **not** a list of people who have consented.

| Commits | Name |
|---|---|
| 128 | Juniormar Organista |
| 69 | Pedro Coimbra |
| 40 | Daniel Garcia |
| 28 | Kainã |
| 4 | Johnatas |
| 3 | Antonio Castelo Filho |

Contributors who worked on the project without committing directly — advisors,
students whose work was merged by someone else, people credited in the papers —
will not show up here. Ask around before assuming this list is complete.
