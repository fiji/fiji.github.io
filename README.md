# The [Fiji Website](https://fiji.sc)

This is the repository for the [Fiji website](https://fiji.sc). It's a simple front-facing website that lists basic information about Fiji, and directs users to more information if they need it. In addition, it has an archive of older bug information, kept for posterity. It is built with [Bootstrap 5](https://getbootstrap.com), [Font Awesome](https://fontawesome.com), and a few other open-source web libraries.

Like Fiji itself, this website is built by the Fiji community, and is continuously improved by contributors all around the world. With that in mind, all forks and PRs are welcome!

## Development Setup

Setting up the website and developing on it is quite simple - all any developer needs is terminal access and a text editor!

First, clone (and fork, if needed) the repository:

```sh
$ git clone https://github.com/fiji/fiji.github.io.git
$ cd fiji.github.io
```

Then, **we have to serve the directory with a HTTP server**. This is because the plugins pane is generated with data from `plugins.json`, which the website gets via an http request. Luckily, there are many ways to start a quick webserver - here is one example with python.

```sh
$ python -m SimpleHTTPServer
Serving HTTP on 0.0.0.0 port 8000 ...
```

And voila, everything should be done! Visit whatever port the HTTP server is opened on (this tends to be `http://localhost:8000/` by default).

## Adding A Plugin

Want to feature a plugin on the Fiji website? Great! All you need to do is submit a PR editing the `plugins.json` file, and optionally uploading a picture.

In `plugins.json`, there is one JSON object titled `plugins`, which is an array of (anonymous) plugin objects. To properly create a plugin object, it needs three mandatory (and one optional) parameters:

* `name`: the name of the plugin. It's alright if it has spaces or other non-URL friendly characters - it's internally slugified when referenced.
* `description`: a brief (one or two line) description of the plugin.
* `link`: a link to the *ImageJ wiki* entry for the plugin, or another informational website if not available.
* `imgUrl`: an (optional) URL to an image that showcases the plugin.
    * The image can either be a relative path (to an image that already exists in this repo, like `site/img/plugins/...`), or an external link to some other resource. However, we prefer that you upload your own image to the repository when possible.
    * The image can be any standard file type that is supported by the web - however, compressed JPGs and PNGs tend to be the best for page performance.
