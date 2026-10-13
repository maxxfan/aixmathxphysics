// View and like counters. GitHub Pages is static, so the counts live on Abacus
// (https://abacus.jasoncameron.dev), a free counter API: /hit increments, /get reads.
(function () {
  var API = "https://abacus.jasoncameron.dev";
  var NS = "maxxfan.github.io";
  var VIEWS = "aixmathxphysics.views";
  var LIKES = "aixmathxphysics.likes";

  var viewsEl = document.getElementById("view-count");
  var likesEl = document.getElementById("like-count");
  var likeBtn = document.getElementById("like-btn");
  if (!viewsEl || !likesEl || !likeBtn) return;

  // Storage can be unavailable (private mode, blocked site data); never let that break the page.
  function read(store, k) {
    try { return window[store].getItem(k); } catch (e) { return null; }
  }
  function write(store, k, v) {
    try { window[store].setItem(k, v); } catch (e) {}
  }

  function fmt(n) {
    return Number(n).toLocaleString("en-US");
  }

  function call(op, key) {
    return fetch(API + "/" + op + "/" + NS + "/" + key).then(function (r) {
      if (r.status === 404 && op === "get") return { value: 0 };
      if (!r.ok) throw new Error(r.status);
      return r.json();
    }).then(function (d) {
      return d.value;
    });
  }

  // Count one view per browser session, so reloading the page doesn't inflate it.
  var seen = read("sessionStorage", "axmp-viewed") === "1";
  call(seen ? "get" : "hit", VIEWS).then(function (v) {
    viewsEl.textContent = fmt(v);
    var label = document.getElementById("view-label");
    if (label) label.textContent = v === 1 ? "view" : "views";
    write("sessionStorage", "axmp-viewed", "1");
  }).catch(function () {
    viewsEl.textContent = "–";
  });

  function setLiked(on) {
    likeBtn.setAttribute("aria-pressed", on ? "true" : "false");
    likeBtn.title = on ? "You liked this" : "Like";
  }

  setLiked(read("localStorage", "axmp-liked") === "1");
  call("get", LIKES).then(function (v) {
    likesEl.textContent = fmt(v);
  }).catch(function () {
    likesEl.textContent = "–";
  });

  // One like per browser: the counter can only go up, so a second click does nothing.
  likeBtn.addEventListener("click", function () {
    if (likeBtn.getAttribute("aria-pressed") === "true") return;
    setLiked(true);
    call("hit", LIKES).then(function (v) {
      likesEl.textContent = fmt(v);
      write("localStorage", "axmp-liked", "1");
    }).catch(function () {
      setLiked(false);
    });
  });
})();
