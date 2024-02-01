# Coding Conventions

This is a list of some things that contributions should adhere to. 
The code still has severe legacy problems with these issues.

1. Keep constructors free of side effects. Prefer using static factory functions.
2. If there are many parameters, use a builder pattern.
3. Never deliberately use `nullptr` values. Prefer `std::optional` instead.
4. Don't use `using namespace std`
5. Use the type system to your advantage.
6. Separate parsing code. In general, learn something about how to separate concerns in programming, learn the SOLID principles, **and** apply them.
7. Use the standard library, including the Standard Template Library. Prefer searching in the C++ standard library over reinventing the wheel.
8. Use the boost library. It is already a dependency. Prefer search in boost over reinventing the wheel.
9. Always try to leave the code in a better (more readable, understandable, maintainable, safer) state than you found it.
10. C++ is hard to read, so don't make it harder than necessary. Code readability is **not optional**. 
    * Use descriptive but concise names for variables, functions, classes, etc.
    * Keep lines short.
    * Prefer vertical lists (e.g. of function arguments) over horizontal lists).
    * Avoid "what" and "how" comments. Prefer "why" comments.
11. If you figure out something really hard and unintuitive, add a comment instead of letting the next programmer figure it out again.
