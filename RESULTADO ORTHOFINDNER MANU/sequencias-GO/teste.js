try {
  var feed = document.getElementById("feed-placeholder");
  console.log("~~Infinite Scroll: Main feed container: ", feed);

  if (feed !== null) {
    var firstObserver = new MutationObserver(function (mutations) {
      mutations.forEach(function (mutation) {
        console.log(
          "~~Primeiro observer: Mutation detectada no #feed-placeholder."
        );

        var bVejaMais = document.querySelector(
          "#feed-placeholder >  div > div > div._l"
        );
        console.log("~~Resultado veja mais: ", bVejaMais);
        if (bVejaMais !== null) {
          console.log("~~Veja mais encontrado! ", bVejaMais);
          firstObserver.disconnect();
          console.log("~~Primeiro observer desconectado");
          var secondObserver = new MutationObserver(function (mutations) {
            mutations.forEach(function (mutation) {
              console.log("~~Setagem do segundo observer");
              window.dataLayer.push({
                event: "pageviewVirtualCS",
                tipoConteudo: "Infinite Scroll",
                urlPathVirtual: document.location.pathname,
              });
            });
          });

          var secondObserverConfig = { childList: true };
          secondObserver.observe(bVejaMais, secondObserverConfig);
        }
      });
    });

    var firstObserverConfig = { childList: true, subtree: true };
    firstObserver.observe(feed, firstObserverConfig);
  } else {
    console.error(
      "~~Erro: Elemento 'feed-placeholder' não encontrado no carregamento inicial."
    );
  }
} catch (e) {
  console.log("~~Erro nos observers:", e);
}

/*try {
  var feed = document.getElementById("feed-placeholder");
  console.log("~~Infinite Scroll: Main feed container: ", feed);

  if (feed !== null) {
    var firstObserver = new MutationObserver(function (mutations) {
      mutations.forEach(function (mutation) {
        console.log("~~Primeiro observer: Mutation detectada no #feed-placeholder.");
        var bVejaMais = document.querySelector("feed-placeholder > div > div > div._l");
        console.log("~~DOM at mutation:", feed);

        if (bVejaMais !== null) {
          console.log("~~Infinite Scroll 2: bVejaMais encontrado: ", bVejaMais);
          firstObserver.disconnect();

          var secondObserver = new MutationObserver(function (mutations) {
            mutations.forEach(function (mutation) {
              console.log("~~Segundo Observer: Mutation detectada no bVejaMais.");
              window.dataLayer.push({
                event: "pageviewVirtualCS",
                tipoConteudo: "Infinite Scroll",
                urlPathVirtual: document.location.pathname,
              });
              console.log("Evento Infinite Scroll disparado!");
            });
          });

          var secondObserverConfig = { childList: true };
          secondObserver.observe(bVejaMais, secondObserverConfig);
          console.log("~~Segundo Observer configurado em bVejaMais.");
        } else {
          console.log("~~Primeiro observer: bVejaMais ainda não encontrado.");
        }
      });
    });

    console.log("~~Primeiro observer: Iniciando observação em 'feed-placeholder'.");
    var firstObserverConfig = { childList: true, subtree: true };
    firstObserver.observe(feed, firstObserverConfig);
  } else {
    console.error("~~Erro: Elemento 'feed-placeholder' não encontrado no carregamento inicial.");
  }
} catch (e) {
  console.error("~~Erro nos observers:", e.message, e.stack);
}*/
