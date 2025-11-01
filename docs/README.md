# Current outline
This is subject to change.
* Introduction
* Background math
    * Linear algebra
        * (This can link to external resources, but readers should be guided on what parts of linear algebra are worth learning, and making these docs self-contained would be a good long-term goal. If we do use external links, we should include a date so that readers know when health of each link was last checked.)
        * Vectors
        * Matrices, matrix-vector multiplication, and its meaning
        * Matrix-matrix multiplication and its meaning
            * (We should likely explain both the "transformation" and "change of basis" interpretations)
        * 3D examples (assuming previous sections have used 2D examples)
        * Dot products
        * Projections, reflections, rotations
        * Homogeneous coordinates and translations
        * (TODO: Check if cross products or determinants are useful for Hypermine)
    * Spherical geometry
        * (We want such a section because it's a good segue to problem solving techniques for hyperbolic geometry problems)
        * Representing points as unit vectors
        * Projections, reflections, rotations, translations (which are rotations expressed differently)
    * Hyperbolic geometry
        * Minkowski space
        * Representing points as "normalized" vectors
        * Projections, reflections, rotations, translations, horo-rotations
        * (Do we put advanced shape-casting math here? Probably not.)

# Design guidelines
* Docs should live in the repository itself rather than a separate wiki.
    * This helps keep docs in sync with the code and allows the quality of docs to be enforced with pull requests. The wiki has a disadvantage of forcing users to make unilateral edits to it.
* Docs should have a suggested linear order to read them.
* There should be a way for readers to know if it's safe to skip a section.
    * This can for instance be done with "After reading this section, readers should be able to (...). This could be questions they can answer.
* Since math is heavily involved, exercises can be useful as knowledge checks.
* While there is a linear order, readers should also be able to tell what parts they can skip if they are reading the docs for a particular purpose.
* We should be prepared to keep placeholders in the docs, potentially with links to external sources, as this would allow us to separate the tasks of writing and organizing documentation.
* Since geometry is heavily involved, the docs should contain pictures.
    * Interactive elements would also be helpful.
* To avoid running up against GitHub limits or making repositories take longer to clone, images and videos will need to be generated on the reader's machine.
    * One option would be to add functionality to Hypermine itself for these visualizations. It is not decided what the preferred approach is.
* We should try to keep the reader interested/motivated, making the documentation enjoyable to read.
    * Animations and interactive visualizations can help with this a lot (in addition to helping the learning process).
